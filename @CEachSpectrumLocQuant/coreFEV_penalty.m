function [abundance,frageffe]=coreFEV_penalty(~,X,massArrangement,penalty_factor)
% Core program for solving the optimization problem, does not depend on individual spectra, variable model of fragmentation efficiency, uses matching score as penalty
% Input:
%   X - the A matrix in $A\beta$
%   massArrangement - the matrix of all combinations of modification masses, each row represents an IMP, each column is the mass shift at each site
%   penalty_factor - the penalty for each IMP
% Output: 
%   abundance - the abundance of each IMP, relative abundance
%   frageffe - fragmentation efficiency of each observed ion
m=size(massArrangement,1);
H=X'*X; % Quadratic term matrix of the optimization function
f=zeros(size(X,2),1);
f(1:m) = penalty_factor;
lb=zeros(size(X,2),1); % Some are >=0, some are >=1
for r=m+1:length(lb)
    lb(r)=1;
end
options = optimoptions('quadprog','Display','off');
[abundance,~,~] = quadprog(H,f,[],[],[],[],lb,[],[],options);
abundance(1:m)=abundance(1:m).*(abs(abundance(1:m))>1e-8); % If less than 1e-8, consider the value as 0
frageffe=1./abundance(m+1:end); % For the latter analysis, fragmentation efficiency of each ion
if abs(sum(abundance(1:m)))>1e-8 % Normalize the obtained quantitative ratio, and avoid normalizing to seemingly reasonable values when all are close to 0
    abundance=abundance(1:m)/sum(abundance(1:m));
else
    abundance=zeros(m,1);
end

% The following is the fragmentation efficiency of each observed ion, which does not affect the discrimination and quantification
frageffe=frageffe.*(abs(frageffe)>1e-8); % If less than 1e-8, consider the value as 0
if abs(sum(frageffe))>1e-8
    frageffe=frageffe/max(frageffe);
else
    frageffe=zeros(length(abundance)-m,1);
end

end