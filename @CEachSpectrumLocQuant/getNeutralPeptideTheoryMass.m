function [theoryMass]=getNeutralPeptideTheoryMass(obj,nonTargetMod)
% Calculate the theoretical mass of a neutral peptide with fixed modifications
% Output: theoryMass is the neutral peptide mass with fixed modifications

% Convert the sequence into numbers using ASCII characters to access the AA mass list, and calculate the mass
vPepAAMass=CConstant.vAAmass(obj.m_pepSeq-64);
theoryMass=sum(vPepAAMass)+2*CConstant.hmass+CConstant.omass;

% Add fixed modifications
for idx=1:size(nonTargetMod,1)
    theoryMass=theoryMass+nonTargetMod{idx,3};
end
end