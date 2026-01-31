function [match_ions] = match_1p_1s(~, peptide, spectrum, tolerance)
% Match one theoretical peptide ions with peaks in the spectrum.
% input:
%   peptide (struct)
%       fields: seq (char), mod_mass (1 x M double), mod_pos (1 x M double)
%       seq:    peptide sequence
%       mod_mass: modification masses
%       mod_pos:  modification positions on the peptide
%   spectrum (struct)
%       fields: peaks (N x 2 double [m/z, intensity]), pre_charge (1 x 1 double/int)
%       peaks:      experimental peaks
%       pre_charge: precursor charge according to the spectrum
%   tolerance (struct)
%       fields: value (double), is_ppm (logical)
%       value:  mass tolerance value
%       is_ppm: true if tolerance is in ppm
% output:
%   match_ions (L x 4 double)
%       [match_type, match_pos, charge, expe_which]
%       match_type: 1 for b-ion, 2 for y-ion
%       match_pos:  ion position on peptide sequence
%       charge:     theoretical ion charge
%       expe_which: index of matched experimental peak

ion_types = [1, 2]; % only consider the b/y ion here
match_ions = zeros(length(ion_types)*length(peptide.seq),4);
idx_match_record = 1;
tol = tolerance.value;

% get theoretical ions of the peptide, [types, positions, charge, m/zs]
theo_ions = get_theoretical_ions(peptide, spectrum.pre_charge);

% match between experiment peaks and theoretical peaks
for iExpPeak = 1:size(spectrum.peaks,1)
    if tolerance.is_ppm
        tol = spectrum.peaks(iExpPeak,1)*tolerance.value / 1e6;
    end
    % find theoretical ions in the tolerance interval
    idx_first_match = find(theo_ions(:,4)>spectrum.peaks(iExpPeak,1)-tol & ...
        theo_ions(:,4)<spectrum.peaks(iExpPeak,1)+tol);
    if ~isempty(idx_first_match)
        % record the first match, its priority is highest
        match_ions(idx_match_record,:) = ...
            [theo_ions(idx_first_match(1),[1,2,3]),iExpPeak];
        idx_match_record = idx_match_record+1;
    end
end
match_ions(idx_match_record:end,:)=[];

% get sorted unique ions
match_ions = sortrows(match_ions,[1,2,3]);
end


function theo_ions = get_theoretical_ions(peptide, pre_charge)
% Get theoretical ions.
% input:
%   peptide (struct)
%       fields: seq (char), mod_mass (1 x M double), mod_pos (1 x M double)
%   pre_charge (1 x 1 double/int)
%       precursor charge according to the spectrum
% output:
%   theo_ions (T x 4 double)
%       theoretical ions: [types, positions, charge, m/zs]

% calculate the max charge of fragment ions according to precursor charge
if pre_charge <= 2
    maxCharge = 1;
else
    maxCharge = 2;
end

% neutral amino acid with modifications
vPepAAMass = [0,CConstant.vAAmass(peptide.seq-64),0]';% Masses of various amino acids without modification
for idx = 1:length(peptide.mod_mass)
    vPepAAMass(peptide.mod_pos(idx)+1) = ...
        vPepAAMass(peptide.mod_pos(idx)+1) + peptide.mod_mass(idx);
end

% Ion_types means: (this order is meaningful)
% 1 y ion, 11 y NH3-loss, 12 y H2O-loss
% 2 b ion, 21 b NH3-loss, 22 b H2O-loss
% 30 a ion, 31 a NH3-loss, 31 a H2O-loss
ion_types = [1, 2]; % only consider the b/y ion here
% ion_types = [1,11,12,2,21,22,30,31,32]; % consider a/b/y ions and their loss
pep_len = length(peptide.seq);
theo_ions = zeros(length(ion_types)*maxCharge*(pep_len-1), 4);
iStart = 1;
for i_charge = 1:maxCharge
    b = (cumsum(vPepAAMass(2:pep_len))+i_charge*CConstant.pmass)/i_charge;
    y = (cumsum(vPepAAMass(pep_len+1:-1:3))+2*CConstant.hmass+CConstant.omass ...
        + i_charge*CConstant.pmass)/i_charge;

    % b ion
    if any(2==ion_types)
        iEnd = iStart+pep_len-2;
        theo_ions(iStart:iEnd,:) = [repmat(2,pep_len-1,1), (1:pep_len-1)',...
            repmat(i_charge,pep_len-1,1), b];
        iStart = iEnd+1;
    end
    % b NH3-loss
    if any(21==ion_types)
        iEnd = iStart+pep_len-2;
        theo_ions(iStart:iEnd,:) = [repmat(21,pep_len-1,1), (1:pep_len-1)',...
            repmat(i_charge,pep_len-1,1), b-(CConstant.nmass+3*CConstant.hmass)/i_charge];
        iStart = iEnd+1;
    end
    % b H2O-loss
    if any(22==ion_types)
        iEnd = iStart+pep_len-2;
        theo_ions(iStart:iEnd,:) = [repmat(22,pep_len-1,1), (1:pep_len-1)',...
            repmat(i_charge,pep_len-1,1), b-(2*CConstant.hmass+CConstant.omass)/i_charge];
        iStart = iEnd+1;
    end

    % y ion
    if any(1==ion_types)
        iEnd = iStart+pep_len-2;
        theo_ions(iStart:iEnd,:) = [ones(pep_len-1,1), (1:pep_len-1)',...
            repmat(i_charge,pep_len-1,1), y];
        iStart = iEnd+1;
    end
    % y NH3-loss
    if any(11==ion_types)
        iEnd = iStart+pep_len-2;
        theo_ions(iStart:iEnd,:) = [repmat(11,pep_len-1,1), (1:pep_len-1)',...
            repmat(i_charge,pep_len-1,1), y-(CConstant.nmass+3*CConstant.hmass)/i_charge];
        iStart = iEnd+1;
    end
    % y H2O-loss
    if any(12==ion_types)
        iEnd = iStart+pep_len-2;
        theo_ions(iStart:iEnd,:) = [repmat(12,pep_len-1,1), (1:pep_len-1)',...
            repmat(i_charge,pep_len-1,1), y-(2*CConstant.hmass+CConstant.omass)/i_charge];
        iStart = iEnd+1;
    end

%     % a ion
%     if any(30==ion_types)
%         iEnd = iStart+pep_len-2;
%         theo_ions(iStart:iEnd,:) = [repmat(3,pep_len-1,1), (1:pep_len-1)',...
%             repmat(i_charge,pep_len-1,1), b-(CConstant.cmass+CConstant.omass)/i_charge];
%         iStart = iEnd+1;
%     end
%     % a NH3-loss
%     if any(31==ion_types)
%         iEnd = iStart+pep_len-2;
%         theo_ions(iStart:iEnd,:) = [repmat(31,pep_len-1,1), (1:pep_len-1)',...
%             repmat(i_charge,pep_len-1,1), b-(CConstant.cmass+CConstant.omass+CConstant.nmass+ ...
%             3*CConstant.hmass)/i_charge];
%         iStart = iEnd+1;
%     end
%     % a H2O-loss
%     if any(32==ion_types)
%         iEnd = iStart+pep_len-2;
%         theo_ions(iStart:iEnd,:) = [repmat(32,pep_len-1,1), (1:pep_len-1)',...
%             repmat(i_charge,pep_len-1,1), b-(CConstant.cmass+CConstant.omass+2*CConstant.hmass+ ...
%             CConstant.omass)/i_charge];
%         iStart = iEnd+1;
%     end
end
theo_ions = sortrows(theo_ions,1); % sort by type, y/b/a/y-loss ...
end