function [abundance,frageff]=coreFEV_lasso(X,massArrangement,lambda)
% coreFEV_lasso - Solve FEV model with lasso regularization.
% Inputs:
%   X (N x P double)
%       Design matrix for FEV model.
%   massArrangement (M x S double)
%      Matrix of all possible mod mass combinations. Each row is a case, each column is the mass shift at each possible modification site.
%   lambda (1 x 1 double)
%       Lasso regularization parameter.
% Outputs:
%   abundance (M x 1 double)
%       Relative abundance for each peptidoform.
%   frageff (Q x 1 double)
%       Normalized fragment efficiency.

m=size(massArrangement,1);
H=X'*X;
vv=ones(size(X,2),1);
f=0.5*lambda*vv;
f(m+1:end)=0;
lb=zeros(size(X,2),1);
for r=m+1:length(lb)
    lb(r)=1;
end
options = optimoptions('quadprog','Display','off');
[abundance,~,~] = quadprog(H,f,[],[],[],[],lb,[],[],options);
abundance(1:m)=abundance(1:m).*(abs(abundance(1:m))>1e-8);
frageff=1./abundance(m+1:end);
if abs(sum(abundance(1:m)))>1e-8
    abundance=abundance(1:m)/sum(abundance(1:m));
else
    abundance=zeros(m,1);
end

% Normalize fragment efficiency coefficients (informational; does not affect abundance discrimination).
frageff=frageff.*(abs(frageff)>1e-8);
if abs(sum(frageff))>1e-8
    frageff=frageff/max(frageff);
else
    frageff=zeros(length(abundance)-m,1);
end

end