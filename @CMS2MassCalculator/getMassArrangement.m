function [ bSuccess,inxSites,massArrangement,warning_msg ] = getMassArrangement(ctx,fixedPosMod)
% Calculate all possible arrangements of multiple modifications and return an ordered mass matrix
% Inputs:
%   ctx (struct)
%       Required fields: m_pepSeq, m_isProtN, m_isProtC, m_variableModNameMass,
%       m_ms1_tolerance, m_dPrecursorMass, m_enzyme, m_strSpecName.
%   fixedPosMod (K x 3 cell)
%       Fixed modification list [position, mod_name, mod_mass].
% Outputs:
%   bSuccess (1 x 1 logical)
%       True when feasible arrangements are found and candidate count is acceptable.
%   inxSites (1 x S double/int)
%       Modification-site positions corresponding to massArrangement columns.
%   massArrangement (M x S double)
%       Candidate peptidoform mass arrangement matrix.
%   warning_msg (1 x 1 char/string)
%       Warning text for no-feasible/too-many-candidates scenarios.

warning_msg = [];
bSuccess = false;

deltamass = ctx.m_dPrecursorMass - CMS2MassCalculator.getNeutralPeptideTheoryMass(ctx,fixedPosMod);

% eachSpecfinVariList is a table of [amino acid, number of modification types with this specificity site, set of positions of this amino acid modification in the user-specified list], for those modifications specified by the user
eachSpecfinVariList=findVariableInSeq(ctx.m_variableModNameMass, ...
    ctx.m_pepSeq,ctx.m_isProtN,ctx.m_isProtC);

% All positions on the sequence where modification can occur
inxSites=[];
% Maximum number of modifications for each specificity
maxNumEachAA=zeros(size(eachSpecfinVariList,1),1);
% Maximum number of each modification
maxNumEachMod=zeros(size(ctx.m_variableModNameMass,1),1);
% The above two must be column vectors

for z=1:size(eachSpecfinVariList,1)
    if isequal(eachSpecfinVariList{z,1},'N-term')
        inxSitesZ=0;
        maxNumEachAA(z)=1;
        maxNumEachMod(eachSpecfinVariList{z,3})=1;
    elseif isequal(eachSpecfinVariList{z,1},'C-term')
        inxSitesZ=length(ctx.m_pepSeq)+1;
        maxNumEachAA(z)=1;
        maxNumEachMod(eachSpecfinVariList{z,3})=1;
    else
        inxSitesZ=strfind(ctx.m_pepSeq,eachSpecfinVariList{z,1});
        maxNumEachAA(z)=length(inxSitesZ);
        maxNumEachMod(eachSpecfinVariList{z,3})=maxNumEachAA(z);
    end
    inxSites=[inxSites,inxSitesZ]; %#ok<AGROW>
end

% Calculate all possible modification combinations, each row is a case, each position represents the number of occurrences of various modifications
if ctx.m_ms1_tolerance.isppm
    ms1_tolerance = ctx.m_ms1_tolerance.value*ctx.m_dPrecursorMass*1e-6;
else
    ms1_tolerance = ctx.m_ms1_tolerance.value;
end

modComb = get_weights_comb([ctx.m_variableModNameMass{:,3}]',maxNumEachMod,deltamass,...
    ms1_tolerance,eachSpecfinVariList,maxNumEachAA);

% Build arrangement matrix from feasible combinations and apply enzyme-specific filtering.
if isempty(modComb) || ~all(sum(modComb,2))
    % Cannot be matched or is an unmodified peptide
    massArrangement=[];
else
    [massArrangement, is_too_many_candidate] = getMassArrangementUsingComb(...
        modComb,ctx.m_variableModNameMass, eachSpecfinVariList,maxNumEachAA);

    if is_too_many_candidate
        massArrangement = [];
        warning_msg = ['There are too many candidate peptidoforms for ',...
            ctx.m_pepSeq, ' in ', ctx.m_strSpecName, '!\n'];
        return
    end

    % The first row is the positions on the sequence, each subsequent row is the mass shift at each position
    sortSitesArrang=(sortrows([inxSites;massArrangement]'))';   % sort by columns
    inxSites=sortSitesArrang(1,:);
    massArrangement=sortSitesArrang(2:end,:);

    % Filtering some unavailable mass arrangement
    % Use enzyme right now
    if isequal(ctx.m_enzyme.name,'trypsin')
        idx_C_term_AA = inxSites==length(ctx.m_pepSeq);
        if any(idx_C_term_AA)
            is_delete_C_term = true(1,size(massArrangement,1));
            is_delete_C_term(massArrangement(:,idx_C_term_AA)<ms1_tolerance) = false;
            for idx_limit = 1:length(ctx.m_enzyme.limits)
                is_delete_C_term((ctx.m_enzyme.limits(idx_limit)-ms1_tolerance<massArrangement(:,idx_C_term_AA)&...
                    massArrangement(:,idx_C_term_AA)<ctx.m_enzyme.limits(idx_limit)+ms1_tolerance)) = false;
            end
           % If a modification occurs in C-term K/R, and the mass is greater
            %   then 14.0266 (average mass of methylation, slightly
            %   greater than the monoisotopic mass of methylation, which
            %   is 14.015650), the PSM is absolutely wrong.
            %         massArrangement(massArrangement(:,idx_C_term_AA)>14.0266,:) = [];
            massArrangement(is_delete_C_term,:) = [];
        end
    end
end

if isempty(massArrangement)
    warning_msg = ['There is no feasible modification configurations for ',...
        ctx.m_pepSeq, ' in ', ctx.m_strSpecName, '!\n'];
elseif size(massArrangement,1) > 5000
    warning_msg = ['There are too many candidate peptidoforms for ',...
        ctx.m_pepSeq, ' in ', ctx.m_strSpecName, '!\n'];
else
    bSuccess = true;
end

end

function all_res = get_weights_comb(each_mass_shift, each_max_num, delta_mass, ...
    tolerance, eachSpecfinVariList, maxNumEachAA)
% Enumerate feasible modification-count combinations with bounded DFS pruning.
% Function for finding all possible combinations of weights in "each_mass_shift" vector,
% such that the sum of selected weights is nearly equal to "delta_mass".
%
% Input arguments:
% - each_mass_shift (R x 1 double): a vector of positive integers representing the available weights.
% - each_max_num (R x 1 double/int): a non-negative integer vector of the same length as "weights", indicating the maximum
%       number of times each weight can be used.
% - delta_mass (1 x 1 double): a positive integer representing the desired total weight.
% - tolerance (1 x 1 double): the error tolerance.
% - eachSpecfinVariList (A x 3 cell): the specificity and the index of the which
%       [each_mass_shift] belongs to this specificity.
% - maxNumEachAA (A x 1 double/int): the total number of the specificity.
% Output arguments:
% - all_res (C x R double): a matrix of all combination result, a row according to one
%       combination.

buff_len = 100;
type_len = length(each_mass_shift);
all_res_idx = 0;
all_res = zeros(buff_len,type_len);

states = repmat(struct('used', zeros(type_len,1), 'index', 1, 'target', delta_mass, 'limits', each_max_num), buff_len, 1);
top = 1;

% Initialize the first state
states(top).target = delta_mass;

while top > 0
    curr_state = states(top);
    top = top - 1;

    % The sum of the each specificity should not be greater than the
    %   limit of each specificity.
    if_skip = false;
    for idx_spec = 1:size(eachSpecfinVariList,1)
        if sum(curr_state.used(eachSpecfinVariList{idx_spec,3})) > maxNumEachAA(idx_spec)
            if_skip = true;
            break;
        end
    end
    if if_skip
        continue;
    end

    % Check if we have found a valid combination
    if abs(curr_state.target) < tolerance
        all_res_idx = all_res_idx + 1;
        if all_res_idx > size(all_res,1)
            all_res(all_res_idx+buff_len,:) = 0;
        end
        all_res(all_res_idx,:) = curr_state.used';
        continue;
    elseif all(curr_state.limits == 0)
        continue;
    elseif ~any(each_mass_shift(curr_state.index:end)<0)...
            && (curr_state.limits(curr_state.index:end)' * each_mass_shift(curr_state.index:end) < ...
            curr_state.target - tolerance || curr_state.target < 0 || ...
            min(each_mass_shift(curr_state.index:end)) > curr_state.target + tolerance)
        % Every available element is greater than 0,
        % so a very big target or a <0 target lead to nothing,
        % and the smallest should not be greater than the target.
        continue;
    elseif ~any(each_mass_shift(curr_state.index:end)>0)...
            && (curr_state.limits(curr_state.index:end)' * each_mass_shift(curr_state.index:end) > ...
            curr_state.target + tolerance || curr_state.target > 0 || ...
            max(each_mass_shift(curr_state.index:end)) < curr_state.target - tolerance)
        % every available element is less than 0,
        % so a very small target or a >0 target or lead to nothing,
        % and the least should not be less than the target.
        continue;
    else
        % Add the next state where the current mass difference is used
        if curr_state.limits(curr_state.index) > 0
            top = top + 1;
            if top > length(states)
                states(top+buff_len) = ...
                    struct('used', zeros(1,type_len), 'index', 1, 'target', delta_mass, 'limits', each_max_num);
            end
            states(top) = curr_state;
            states(top).limits(curr_state.index) = states(top).limits(curr_state.index) - 1;
            states(top).used(curr_state.index) = states(top).used(curr_state.index) + 1;
            states(top).target = curr_state.target - each_mass_shift(curr_state.index);
        end

        % Add the next state where the current mass difference is not used
        if curr_state.index < length(each_mass_shift)
            top = top + 1;
            if top > length(states)
                states(top+buff_len) = ...
                    struct('used', zeros(1,type_len), 'index', 1, 'target', delta_mass, 'limits', each_max_num);
            end
            states(top) = curr_state;
            states(top).index = states(top).index + 1;
        end
    end
end
all_res(all_res_idx+1:end,:) = [];
end

function [massArrangement, is_too_many_candidate]=getMassArrangementUsingComb(modComb,variModNameMass, ...
    eachSpecfinVariList,maxNumEachAA)
% Expand each combination into site-level permutations of mass assignments.
% Calculate all combinations of modification mass + modification sites (just the order of potential modifications on the sequence, not the actual positions)
% Inputs:
%   modComb (C x R double) - All possible modification combinations. Each row is one combination, the column order is consistent with the user-specified modification list
%   variModNameMass (R x 3 cell) - a matrix of modification types and modification masses, consistent with the order of all modifications specified by the user
%   eachSpecfinVariList (A x 3 cell) - a table of [amino acid, number of modification types with this specificity site, positions of this amino acid modification in the user-specified list]
%   maxNumEachAA (A x 1 double/int) - the number of positions where various amino acids may be modified on the peptide sequence
% Outputs:
%   massArrangement (M x S double) - a matrix of various combinations of modification masses, each row is a case, each column is the mass shift at several possible modification sites, the columns are organized by amino acids (block matrix) and cannot be used directly, some processing is needed.
%   is_too_many_candidate (1 x 1 logical) - whether there are too many candidate peptidoforms

massArrangement=[];
is_too_many_candidate = false;
for idxComb=1:size(modComb,1)
    massArraEachComb=[];
    for idxAA=1:size(eachSpecfinVariList,1)
        maxNumCurAA=maxNumEachAA(idxAA);
        modPosesVariList=eachSpecfinVariList{idxAA,3};
        if maxNumCurAA==0
            continue;
        else
            massConfigH=zeros(1,maxNumCurAA);
            % The indices of these existing modifications, in the order of the user-specified modification list
            inx=find(modComb(idxComb,modPosesVariList)~=0);
            tmp=0;
            for j=1:length(inx)
                if modComb(idxComb,modPosesVariList(inx(j)))>1
                    for r=1:modComb(idxComb,modPosesVariList(inx(j)))
                        tmp=tmp+1;
                        massConfigH(1,tmp)=variModNameMass{modPosesVariList(inx(j)),3};
                    end
                elseif modComb(idxComb,modPosesVariList(inx(j)))==1
                    tmp=tmp+1;
                    massConfigH(1,tmp)=variModNameMass{modPosesVariList(inx(j)),3};
                end
            end
        end

        if length(massConfigH)>10
            is_too_many_candidate = true;
            return
        end

        massArraEachAA=perms(massConfigH);
        massArraEachAA=unique(massArraEachAA,'rows');

        if isempty(massArraEachComb)
            massArraEachComb=massArraEachAA;
        else
            tmpMassPailie=[];
            for f=1:size(massArraEachComb,1)
                for g=1:size(massArraEachAA,1)
                    tmpMassPailie0=[massArraEachComb(f,:),massArraEachAA(g,:)];
                    tmpMassPailie=[tmpMassPailie;tmpMassPailie0]; %#ok<AGROW>
                end
            end
            massArraEachComb=tmpMassPailie;
        end
    end
    massArrangement=[massArrangement;massArraEachComb]; %#ok<AGROW>
end
end

function [eachSpecfinVariList]=findVariableInSeq(variableModNameMass,pepSeq,isProtN,isProtC)
% Group variable modifications by specificity (N-term, A..Z, C-term).
% Extract various modification types from a long string of input modifications
% Inputs:
%    variableModNameMass - the modification type
%    pepSeq - the peptide sequence
%    isProtN - whether the peptide sequence is a protein N-terminal
%    isProtC - whether the peptide sequence is a protein C-terminal
% Outputs:
%    eachSpecfinVariList - a table of various amino acids and modifications, specifically a 3-column matrix. Each row is an amino acid, 
%        the first column is the amino acid type,
%        the second column is the number of modification types for the same amino acid, 
%        the third column is the set of positions of the same amino acid modification in the user-specified modification string

% Note: Cannot handle the case where a specified modification contains multiple amino acids

% Establish a list in the order of N-term,A,B,C,D,E,F,G,...,X,Y,Z,C-term, and delete the unused ones at the end
eachSpecfinVariList=cell(28,3);
eachSpecfinVariList(:,2) = {0};
for idx=1:size(variableModNameMass,1)
    strSpecfct=variableModNameMass{idx,2};

    % Consider N/C-terminus and ordinary amino acids separately, not considering the case where there are multiple amino acids in one
    if contains(strSpecfct,'N-term')
        cellTemp=split(strSpecfct,'N-term');
        if isequal(cellTemp{1},'Protein') && ~isProtN
            continue;
        end
        if ~isempty(cellTemp{2})
            if ~startsWith(pepSeq,cellTemp{2})
                % N-terminal specificity, and there is an amino acid restriction, skip if not satisfied
                continue;
            end
        end
        eachSpecfinVariList{1,1} = 'N-term';
        eachSpecfinVariList{1,2} = eachSpecfinVariList{1,2}+1;
        eachSpecfinVariList{1,3} = [eachSpecfinVariList{1,3},idx];
    elseif contains(strSpecfct,'C-term')
        cellTemp=split(strSpecfct,'C-term');
        if isequal(cellTemp{1},'Protein') && ~isProtC
            continue;
        end
        if ~isempty(cellTemp{2})
            if ~endsWith(pepSeq,cellTemp{2})
                % C-terminal specificity, and there is an amino acid restriction, skip if not satisfied
                continue;
            end
        end
        eachSpecfinVariList{28,1} = 'C-term';
        eachSpecfinVariList{28,2} = eachSpecfinVariList{28,2}+1;
        eachSpecfinVariList{28,3} = [eachSpecfinVariList{28,3},idx];
    else
        idxList=strSpecfct-'A'+2;
        eachSpecfinVariList{idxList,1} = strSpecfct;
        eachSpecfinVariList{idxList,2} = eachSpecfinVariList{idxList,2}+1;
        eachSpecfinVariList{idxList,3} = [eachSpecfinVariList{idxList,3},idx];
    end
end
eachSpecfinVariList(cellfun(@isempty,eachSpecfinVariList(:,1)),:)=[];
end