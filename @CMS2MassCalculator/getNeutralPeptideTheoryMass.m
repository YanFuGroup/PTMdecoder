function [theoryMass]=getNeutralPeptideTheoryMass(ctx,nonTargetMod)
% Calculate the theoretical mass of a neutral peptide with fixed modifications
% Inputs:
%   ctx (struct)
%       Required field: m_pepSeq.
%   nonTargetMod (K x 3 cell)
%       Fixed modification list [position, mod_name, mod_mass].
% Outputs:
%   theoryMass (1 x 1 double)
%       Neutral peptide theoretical mass including fixed modifications.

vPepAAMass=CConstant.vAAmass(ctx.m_pepSeq-64);
theoryMass=sum(vPepAAMass)+2*CConstant.hmass+CConstant.omass;

% Add fixed modifications
for idx=1:size(nonTargetMod,1)
	theoryMass=theoryMass+nonTargetMod{idx,3};
end
end