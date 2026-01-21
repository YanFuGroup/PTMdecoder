function [ bSuccess,inxSites,massArrangement,warning_msg ] = getMassArrangement(obj,fixedPosMod)
% Calculate all possible arrangements of multiple modifications, represented as mass, and return a ordered mass matrix
% Input: 
%   fixedPosMod - the fixed modification, which is added before calculating the modification combination
% Output: 
%   bSuccess - whether successful
%   inxSites - the actual positions of modifications on the peptide for each column of massArrangement
%   massArrangement - the various arrangements of modifications, given as a matrix in the form of mass.
%   Each row - a possible modification arrangement, each column represents a position
%   warning_msg - the warning message

warning_msg = [];
bSuccess = false;

% Calculate the mass shift contributed by modifications, experimental minus theoretical (with variable modifications) precursor mass
deltamass=obj.m_dPrecursorMass-obj.getNeutralPeptideTheoryMass(fixedPosMod);


% eachSpecfinVariList is a table of [amino acid, number of modification types with this specificity site, set of positions of this amino acid modification in the user-specified list], for those modifications specified by the user
eachSpecfinVariList=findVariableInSeq(obj.m_variableModNameMass, ...
    obj.m_pepSeq,obj.m_isProtN,obj.m_isProtC);

inxSites=[];  % All positions on the sequence where modification can occur
maxNumEachAA=zeros(size(eachSpecfinVariList,1),1);       % Maximum number of modifications for each specificity
maxNumEachMod=zeros(size(obj.m_variableModNameMass,1),1);% Maximum number of each modification
% The above two must be column vectors

for z=1:size(eachSpecfinVariList,1)
    if isequal(eachSpecfinVariList{z,1},'N-term')
        inxSitesZ=0;
        maxNumEachAA(z)=1;
        maxNumEachMod(eachSpecfinVariList{z,3})=1;
    elseif isequal(eachSpecfinVariList{z,1},'C-term')
        inxSitesZ=length(obj.m_pepSeq)+1;
        maxNumEachAA(z)=1;
        maxNumEachMod(eachSpecfinVariList{z,3})=1;
    else
        inxSitesZ=strfind(obj.m_pepSeq,eachSpecfinVariList{z,1});
        maxNumEachAA(z)=length(inxSitesZ);
        maxNumEachMod(eachSpecfinVariList{z,3})=maxNumEachAA(z);
    end
    inxSites=[inxSites,inxSitesZ];% Summarize all possible modification sites
end

% Calculate all possible modification combinations, each row is a case, each position represents the number of occurrences of various modifications
if obj.m_ms1_tolerance.isppm
    ms1_tolerance = obj.m_ms1_tolerance.value*obj.m_dPrecursorMass*1e-6;
else
    ms1_tolerance = obj.m_ms1_tolerance.value;
end
% modComb=getModComb_multi([obj.m_variableModNameMass{:,3}]',deltmass,eachSpecfinVariList, ...
%     maxNumEachAA,maxNumEachMod,ms1_tolerance);
modComb = get_weights_comb([obj.m_variableModNameMass{:,3}]',maxNumEachMod,deltamass,...
    ms1_tolerance,eachSpecfinVariList,maxNumEachAA);


if isempty(modComb) || ~all(sum(modComb,2)) % Cannot be matched or is an unmodified peptide
    massArrangement=[];
else
    % If there are any possible modification combinatioin, use it to establish the modification arrangement
    [massArrangement, is_too_many_candidate] = getMassArrangementUsingComb(...
        modComb,obj.m_variableModNameMass, eachSpecfinVariList,maxNumEachAA);

    if is_too_many_candidate
        massArrangement = [];
        warning_msg = ['There are too many candidate peptidoforms for ',...
            obj.m_pepSeq, ' in ', obj.m_strSpecName, '!\n'];
        return
    end

    % The first row is the positions on the sequence, each subsequent row is the mass shift at each position
    sortSitesArrang=(sortrows([inxSites;massArrangement]'))';% Sort the above matrix by columns
    inxSites=sortSitesArrang(1,:);
    massArrangement=sortSitesArrang(2:end,:);

    % Filtering some unavailable mass arrangement
    % Use enzyme right now
    if isequal(obj.m_enzyme.name,'trypsin')
        % Only these three
        idx_C_term_AA = inxSites==length(obj.m_pepSeq);
        if any(idx_C_term_AA)
            is_delete_C_term = true(1,size(massArrangement,1)); % 1*size(massArrangement,1) true
            is_delete_C_term(massArrangement(:,idx_C_term_AA)<ms1_tolerance) = false;
            for idx_limit = 1:length(obj.m_enzyme.limits)
                is_delete_C_term((obj.m_enzyme.limits(idx_limit)-ms1_tolerance<massArrangement(:,idx_C_term_AA)&...
                    massArrangement(:,idx_C_term_AA)<obj.m_enzyme.limits(idx_limit)+ms1_tolerance)) = false;
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
        obj.m_pepSeq, ' in ', obj.m_strSpecName, '!\n'];
elseif size(massArrangement,1) > 5000
    warning_msg = ['There are too many candidate peptidoforms for ',...
        obj.m_pepSeq, ' in ', obj.m_strSpecName, '!\n'];
else
    bSuccess = true;
end

end

function [modComb]=getModComb_multi(massShiftEachMod,deltmass,eachSpecfinVariList,maxNumEachAA,maxNumEachMod,ms1_tolerance)
% Use deltamass to calculate all possible modification combinations
% Input: 
%   massShiftEachMod - the mass of the user-specified modifications, must be a column vector
%   deltmass - precursor delta mass (experimental minus theoretical)
%   eachSpecfinVariList - a table of [amino acid, number of modification types with this specificity site, positions of this amino acid in the user-specified list]
%   maxNumEachAA - the number of positions where various amino acids may be modified
%   maxNumEachMod - the number of positions where various modifications may occur, all modifications are the maximum number
%   ms1_tolerance - precursor matching error
% Output: 
%   modComb - all possible modification combinations that explain deltamass

% Get a combination matrix of possible numbers of various modifications
modComb = initModComb_multi(maxNumEachMod);

% Filter using distance
modComb = modCombFilter_dist(modComb,massShiftEachMod,deltmass,ms1_tolerance);

% Filter using the number of modifications
modComb = modCombFilter_modNum(modComb,eachSpecfinVariList,maxNumEachAA);
end

function all_res = get_weights_comb(each_mass_shift, each_max_num, delta_mass, ...
    tolerance, eachSpecfinVariList, maxNumEachAA)
% Function for finding all possible combinations of weights in "each_mass_shift" vector,
% such that the sum of selected weights is nearly equal to "delta_mass".
%
% Input arguments:
% - each_mass_shift: a vector of positive integers representing the available weights.
% - each_max_num: a non-negative integer vector of the same length as "weights", indicating the maximum
%       number of times each weight can be used.
% - delta_mass: a positive integer representing the desired total weight.
% - tolerance: the error tolerance.
% - eachSpecfinVariList: the specificity and the index of the which
%       [each_mass_shift] belongs to this specificity.
% - maxNumEachAA: the total number of the specificity.
% Output arguments:
% - all_res: a matrix of all combination result, a row according to one
%       combination.

% Initial
buff_len = 100;
type_len = length(each_mass_shift);
all_res_idx = 0;
all_res = zeros(buff_len,type_len);

% Create an array to store the states
states = repmat(struct('used', zeros(type_len,1), 'index', 1, 'target', delta_mass, 'limits', each_max_num), buff_len, 1);
top = 1;

% Initialize the first state
states(top).target = delta_mass;

% Loop over the states
while top > 0
    % Get the current state from the top of the stack
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

        % Check if we need to skip this state, check the pruning conditions
    elseif all(curr_state.limits == 0)
        % no available mass shift
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

        % Add new state
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
% Return all valid combinations found
end

function modComb=initModComb_multi(maxNumEachMod)
% Use the maximum possible number of modifications occurrences to establish a combination matrix of possible numbers
% Input: 
%   maxNumEachMod - the number of positions where various modifications may occur
% Output: 
%   modComb - a combination matrix of possible numbers of various modifications for subsequent attempts
if length(maxNumEachMod)==1
    modComb=(0:maxNumEachMod)';
else
    p=0:maxNumEachMod(1);
    q=0:maxNumEachMod(2);
    [x,y]=meshgrid(p,q);
    % Establish combinations and merge
    modComb=[x(:),y(:)];
    if length(maxNumEachMod)>=3
        % When there are more than 3, append an increasing number to the end of each previous matrix, recursively increasing the dimension
        for l=3:length(maxNumEachMod)
            modCombTemp=[];
            for i=0:maxNumEachMod(l)
                modCombTemp=[modCombTemp;modComb,repmat(i,size(modComb,1),1)];
            end
            modComb=modCombTemp;
        end
    end
end
end

function [ modComb ] = modCombFilter_dist( modComb,massShiftEachMod,deltmass,ms1_tolerance )
% Filter out modification combinations that do not meet the criteria using distance
% Input: 
%   modComb - a combination matrix of possible numbers of various modifications
%   massShiftEachMod - the mass of various modifications, must be a column vector
%   deltmass - precursor mass shift, experimental minus theoretical
%   ms1_tolerance - precursor ion matching error
% Output: 
%   modComb - a combination matrix of possible numbers of various modifications

% Only keep columns with distances less than the threshold
modComb = modComb(abs(modComb*massShiftEachMod-deltmass) <= ms1_tolerance,:);

end

function [ modComb ] = modCombFilter_modNum( modComb,eachSpecfinVariList,maxNumEachAA )
% Filter possible modification combinations using the number of modifications
% Input: 
%   modComb - a combination matrix of possible numbers of various modifications
%   eachSpecfinVariList - a table of [amino acid, number of modification types with this specificity site, positions of this amino acid modification in the user-specified list]
%   maxNumEachAA - the number of positions where various amino acids may be modified, a column vector
% Output: 
%   modComb - a combination matrix of possible numbers of various modifications

% Filter once for each target amino acid modification
for idxModAA=1:size(eachSpecfinVariList,1)
    % Sum the positions corresponding to the amino acid in a row, which cannot exceed the maximum possible number of modifications on the amino acid
    modComb = modComb(sum(modComb(:,eachSpecfinVariList{idxModAA,3}),2)<=maxNumEachAA(idxModAA),:);
end

end


function [massArrangement, is_too_many_candidate]=getMassArrangementUsingComb(modComb,variModNameMass, ...
    eachSpecfinVariList,maxNumEachAA)
% Calculate all combinations of modification mass + modification sites (just the order of potential modifications on the sequence, not the actual positions)
% Input: 
%   modComb - All possible modification combinations. Each row is one combination, the column order is consistent with the user-specified modification list
%   variModNameMass - a matrix of modification types and modification masses, consistent with the order of all modifications specified by the user
%   eachSpecfinVariList - a table of [amino acid, number of modification types with this specificity site, positions of this amino acid modification in the user-specified list]
%   maxNumEachAA - the number of positions where various amino acids may be modified on the peptide sequence
% Output: 
%   massArrangement - a matrix of various combinations of modification masses, each row is a case, each column is the mass shift at several possible modification sites, the columns are organized by amino acids (block matrix) and cannot be used directly, some processing is needed.
%   is_too_many_candidate - whether there are too many candidate peptidoforms
massArrangement=[];
is_too_many_candidate = false;
for idxComb=1:size(modComb,1)
    massArraEachComb=[];
    for idxAA=1:size(eachSpecfinVariList,1) % Traverse various specificities

        maxNumCurAA=maxNumEachAA(idxAA);% The number of positions where this amino acid may be modified on the peptide sequence
        modPosesVariList=eachSpecfinVariList{idxAA,3};
        if maxNumCurAA==0
            continue;
        else
            massConfigH=zeros(1,maxNumCurAA);
            inx=find(modComb(idxComb,modPosesVariList)~=0);% The indices of these existing modifications, in the order of the user-specified modification list
            tmp=0;
            for j=1:length(inx)
                if modComb(idxComb,modPosesVariList(inx(j)))>1
                    for r=1:modComb(idxComb,modPosesVariList(inx(j)))
                        tmp=tmp+1;
                        massConfigH(1,tmp)=variModNameMass{modPosesVariList(inx(j)),3};% Record the mass of each modification
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
        
        massArraEachAA=perms(massConfigH);% Use combinations to get all permutations
        massArraEachAA=unique(massArraEachAA,'rows');% Remove duplicates by row


        if isempty(massArraEachComb)
            massArraEachComb=massArraEachAA;
        else
            % Multiple amino acids can be modified
            tmpMassPailie=[];
            for f=1:size(massArraEachComb) % by row
                for g=1:size(massArraEachAA)
                    tmpMassPailie0=[massArraEachComb(f,:),massArraEachAA(g,:)];
                    tmpMassPailie=[tmpMassPailie;tmpMassPailie0];
                end
            end
            massArraEachComb=tmpMassPailie;% Record with massArraEachComb for the next combination into a high-dimensional modification matrix
        end
    end
    % All combinations of modification types + modification sites
    massArrangement=[massArrangement;massArraEachComb];
end
end

function [eachSpecfinVariList]=findVariableInSeq(variableModNameMass,pepSeq,isProtN,isProtC)
% Extract various modification types from a long string of input modifications
% Input: 
%    variableModNameMass - the modification type
%    pepSeq - the peptide sequence
%    isProtN - whether the peptide sequence is a protein N-terminal
%    isProtC - whether the peptide sequence is a protein C-terminal
% Output: 
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