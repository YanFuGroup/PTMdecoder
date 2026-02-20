function [abundance]=coreGuan_penalty(X,Y,penalty_factor)
% coreGuan_penalty - Guan model solved with peptidoform-dependent penalty.
% Inputs:
%   X (N x P double)
%       X matrix in $Y=X\alpha+\epsilon$
%   Y (N x 1 double)
%       Y vector in $Y=X\alpha+\epsilon$
%   penalty_factor (P x 1 double)
%       penalty factor
% Outputs:
%   abundance (P x 1 double)
%       \alpha vector in $Y=X\alpha+\epsilon$ in the paper, representing the relative abundance of each IMP
H=X'*X;
f=-X'*Y+penalty_factor;
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
