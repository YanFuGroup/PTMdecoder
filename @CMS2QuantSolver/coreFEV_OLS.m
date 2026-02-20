function [abundance]=coreFEV_OLS(X,massArrangement)
% coreFEV_OLS - FEV model solved by constrained OLS (quadprog form).
% Inputs:
%   X (N x P double)
%       The A matrix in $A\beta$, where the first M columns correspond to the M peptidoforms, and the remaining columns correspond to the fragment-efficiency variables.
%   massArrangement (M x S double)
%       Matrix of all possible mod mass combinations. Each row is a case, each column is the mass shift at each possible modification site.
% Outputs:
%   abundance (M x 1 double)
%       Normalized IMP abundance.

H=X'*X;
m=size(massArrangement,1);
lb=zeros(size(X,2),1);
for r=m+1:length(lb)
	lb(r)=1;
end
options = optimoptions('quadprog','Display','off');
[alpha,~,~] = quadprog(H,[],[],[],[],[],lb,[],[],options);
abundance=alpha(1:m).*(abs(alpha(1:m))>1e-8);
if abs(sum(abundance))>1e-8
	abundance = abundance/sum(abundance);
else
	abundance = zeros(m,1);
end
end