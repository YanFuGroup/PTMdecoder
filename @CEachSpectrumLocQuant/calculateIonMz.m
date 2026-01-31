function theoryMz=calculateIonMz(obj,fixedPosMod)
% Calculate the m/z of b and y ions with only fixed modifications
% Input: 
%   fixedPosMod (M x 3 cell)
%       fixed modification list [position, name, mass]
% Output: 
%   theoryMz (L x C double) m/z
%       m/z of various b and y ions

% Maximum charge of fragment ions, prescribed rules
if obj.m_iCharge<=2
    maxCharge=1;
elseif obj.m_iCharge==3
    maxCharge=2;
else
    maxCharge=3;
end

vPepAAMass=[0,CConstant.vAAmass(obj.m_pepSeq-64),0];% Mass of amino acids without modification
for idx=1:size(fixedPosMod,1)
    vPepAAMass(fixedPosMod{idx,1}+1) = vPepAAMass(fixedPosMod{idx,1}+1) + fixedPosMod{idx,3};
end

iPeplen=length(obj.m_pepSeq);
theoryMz=zeros(2*maxCharge,(iPeplen-1));% The number of b and y ions is one less than the number of amino acids in the peptide sequence
for i=1:maxCharge
    b=(cumsum(vPepAAMass(1:iPeplen))+i*CConstant.pmass)/i;
    y=(cumsum(vPepAAMass(iPeplen+2:-1:3))+2*CConstant.hmass+CConstant.omass+i*CConstant.pmass)/i;
    theoryMz(i,:)=b(2:end);% Rows 1 to charge are m/z of b ions
    theoryMz((i+maxCharge),:)=y(2:end);% Rows maxCharge+1 to 2*charge are m/z of y ions
end
theoryMz=theoryMz';% Transposed, each row is a position (b_i and y_i), and each column is a charge state
end