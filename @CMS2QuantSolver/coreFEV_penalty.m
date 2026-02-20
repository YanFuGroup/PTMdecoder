function [abundance,frageffe]=coreFEV_penalty(X,massArrangement,penalty_factor)
% Core solver for variable-fragmentation model using penalty factors
% Input:
%   X (N x P double)
%       Design matrix for variable fragmentation efficiency model.
%   massArrangement (M x S double)
%      Matrix of all possible mod mass combinations. Each row is a case, each column is the mass shift at each possible modification site.
%   penalty_factor (M x 1 double)
%       Penalty factor for each peptidoform, derived from matching scores.
% Output:
%   abundance (M x 1 double)
%       Relative abundance for each peptidoform.
%   frageffe (Q x 1 double)
%       Normalized fragment efficiency.

m=size(massArrangement,1);
H=X'*X;
f=zeros(size(X,2),1);
f(1:m) = penalty_factor;
lb=zeros(size(X,2),1);
for r=m+1:length(lb)
    lb(r)=1;
end
options = optimoptions('quadprog','Display','off');
[abundance,~,~] = quadprog(H,f,[],[],[],[],lb,[],[],options);
abundance(1:m)=abundance(1:m).*(abs(abundance(1:m))>1e-8);
frageffe=1./abundance(m+1:end);
if abs(sum(abundance(1:m)))>1e-8
    abundance=abundance(1:m)/sum(abundance(1:m));
else
    abundance=zeros(m,1);
end

% Below is the fragmentation efficiency of each ion, which does not affect the discrimination and quantification of the model
frageffe=frageffe.*(abs(frageffe)>1e-8);
if abs(sum(frageffe))>1e-8
    frageffe=frageffe/max(frageffe);
else
    frageffe=zeros(length(abundance)-m,1);
end

end