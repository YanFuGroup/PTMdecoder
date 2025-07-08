function lfMass = getPeptideMass(obj, pep_seq)
% Get the mass of each peptide
% Input:
%   pep_seq - The peptide sequence
% Output:
%   lfMass - The mass of the peptide
    
    % Add the mass of each amino acid
    lfMass = sum(CConstant.vAAmass(pep_seq-'A'+1));
    
    % Add the mass of water
    lfMass = lfMass + CConstant.hmass*2 + CConstant.omass;
end