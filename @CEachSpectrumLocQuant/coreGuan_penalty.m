function [abundance]=coreGuan_penalty(~,X,Y,penalty_factor)
% Core program for solving the optimization problem, does not depend on individual spectra, constant model of fragmentation efficiency , uses identification score as penalty
% Input: 
%   X - the X matrix in $Y=X\alpha+\epsilon$
%   Y - the Y vector in $Y=X\alpha+\epsilon$
%   penalty_factor - the penalty factor
% Output: 
%   abundance - the \alpha vector in $Y=X\alpha+\epsilon$ in the paper, representing the relative abundance of each IMP
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
