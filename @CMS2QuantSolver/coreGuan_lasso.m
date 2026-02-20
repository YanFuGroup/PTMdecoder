function [abundance]=coreGuan_lasso(X,Y,lambda)
% coreGuan_lasso - Guan model solved with lasso-style linear penalty.
% Inputs:
%   X (N x P double)
%       X matrix in $Y=X\alpha+\epsilon$
%   Y (N x 1 double)
%       Y vector in $Y=X\alpha+\epsilon$
%   lambda (1 x 1 double)
%       lasso penalty factor
% Outputs:
%   abundance (P x 1 double)
%       \alpha vector in $Y=X\alpha+\epsilon$, representing the relative abundance of each IMP
H=X'*X;
vv=ones(size(X,2),1);
f=-X'*Y+0.5*lambda*vv;
lb=zeros(size(X,2),1);
options = optimoptions('quadprog','Display','off');
[beta,~,~] = quadprog(H,f,[],[],[],[],lb,[],[],options);
beta=beta.*(abs(beta)>1e-8);
if abs(sum(beta))>1e-8
	abundance=beta/sum(beta);
else
	abundance = zeros(length(beta),1);
end
end
