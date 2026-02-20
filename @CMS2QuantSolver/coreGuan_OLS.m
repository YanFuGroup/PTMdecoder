function [abundance]=coreGuan_OLS(X,Y)
% coreGuan_OLS - Guan model solved by constrained OLS.
% Inputs:
%   X (N x P double)
%       X matrix in $Y=X\alpha+\epsilon$
%   Y (N x 1 double)
%       Y vector in $Y=X\alpha+\epsilon$
% Outputs:
%   abundance (P x 1 double)
%       \alpha vector in $Y=X\alpha+\epsilon$, representing the relative abundance of each IMP
H=X'*X;
f=-X'*Y;
lb=zeros(size(X,2),1);
options = optimoptions('quadprog','Display','off');
[abundance,~,~] = quadprog(H,f,[],[],[],[],lb,[],[],options);
abundance = abundance.*(abs(abundance)>1e-8);
if abs(sum(abundance))>1e-8
	abundance = abundance/sum(abundance);
else
	abundance = zeros(length(abundance),1);
end
end