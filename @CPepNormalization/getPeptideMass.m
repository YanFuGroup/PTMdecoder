function lfMass = getPeptideMass(~, pep_seq)
% Get the mass of each peptide
% Input:
%   pep_seq (1 x 1 char/string)
%       Peptide sequence
% Output:
%   lfMass (1 x 1 double) Da
%       Mass of the peptide
    
    % Add the mass of each amino acid
    lfMass = sum(CConstant.vAAmass(pep_seq-'A'+1));
    
    % Add the mass of water
    lfMass = lfMass + CConstant.hmass*2 + CConstant.omass;
end