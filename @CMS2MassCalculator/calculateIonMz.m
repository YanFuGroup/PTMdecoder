function theoryMz=calculateIonMz(ctx,fixedPosMod)
% Calculate the m/z of b and y ions with only fixed modifications
% Input:
%   ctx (struct)
%       Required fields: m_pepSeq, m_iCharge.
%   fixedPosMod (K x 3 cell)
%       Fixed modification list [position, mod_name, mod_mass].
% Output:
%   theoryMz (L x C double)
%       Theoretical fragment m/z matrix. Rows are positions, columns are b/y charge channels.

if ctx.m_iCharge<=2
	maxCharge=1;
elseif ctx.m_iCharge==3
	maxCharge=2;
else
	maxCharge=3;
end

vPepAAMass=[0,CConstant.vAAmass(ctx.m_pepSeq-64),0];% Mass of amino acids without modification
for idx=1:size(fixedPosMod,1)
	vPepAAMass(fixedPosMod{idx,1}+1) = vPepAAMass(fixedPosMod{idx,1}+1) + fixedPosMod{idx,3};
end

pep_length=length(ctx.m_pepSeq);
theoryMz=zeros(2*maxCharge,(pep_length-1));% The number of b and y ions is one less than the number of amino acids in the peptide sequence
for i=1:maxCharge
	b=(cumsum(vPepAAMass(1:pep_length))+i*CConstant.pmass)/i;
	y=(cumsum(vPepAAMass(pep_length+2:-1:3))+2*CConstant.hmass+CConstant.omass+i*CConstant.pmass)/i;
	theoryMz(i,:)=b(2:end);% Rows 1 to charge are m/z of b ions
	theoryMz((i+maxCharge),:)=y(2:end);% Rows maxCharge+1 to 2*charge are m/z of y ions
end
theoryMz=theoryMz';% Transposed, each row is a position (b_i and y_i), and each column is a charge state
end
