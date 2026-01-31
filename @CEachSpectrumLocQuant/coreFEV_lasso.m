function [abundance,frageffe]=coreFEV_lasso(~,X,massArrangement,lambda)
% Core program for solving the optimization problem, does not depend on individual spectra, variable model of fragmentation efficiency, solves using lasso
% Input:
%   X (N x P double) - the A matrix in $A\beta$
%   massArrangement (M x S double) - all possible mod mass combinations matrix. Each row is a case, each column is the mass shift at each possible modification site
%   lambda (1 x 1 double) - lasso regularization parameter
% Output: 
%   abundance (M x 1 double) - $\beta$ in $A\beta$
%   frageffe (P-M x 1 double) - fragmentation efficiency of each ion
m=size(massArrangement,1);
H=X'*X;
vv=ones(size(X,2),1);
f=0.5*lambda*vv;
f(m+1:end)=0;
lb=zeros(size(X,2),1); % some >= 0, some >= 1
for r=m+1:length(lb)
    lb(r)=1;
end
options = optimoptions('quadprog','Display','off');
[abundance,~,~] = quadprog(H,f,[],[],[],[],lb,[],[],options);
abundance(1:m)=abundance(1:m).*(abs(abundance(1:m))>1e-8);  % If less than 1e-8, the value is considered 0
frageffe=1./abundance(m+1:end); % Fragmentation efficiency of each ion
if abs(sum(abundance(1:m)))>1e-8
    abundance=abundance(1:m)/sum(abundance(1:m));
else
    abundance=zeros(m,1);
end

% Below is the fragmentation efficiency of each ion, does not affect the discrimination and quantification of the model
frageffe=frageffe.*(abs(frageffe)>1e-8);    % If less than 1e-8, the value is considered 0
if abs(sum(frageffe))>1e-8
    frageffe=frageffe/max(frageffe);
else
    frageffe=zeros(length(abundance)-m,1);
end

end