function [abundance]=coreFEV_OLS(~,X,massArrangement)
% Core program for solving the optimization problem, does not depend on individual spectra, variable model of fragmentation efficiency, solves using OLS
% Input:
%   X - the A matrix in $A\beta$
%   massArrangement - all possible mod mass combinations matrix. Each row is a case, each column is the mass shift at each possible modification site
% Output: 
%   abundance - $\beta$ in $A\beta$
H=X'*X;
m=size(massArrangement,1);
lb=zeros(size(X,2),1);
for r=m+1:length(lb)
    lb(r)=1;
end
options = optimoptions('quadprog','Display','off');
[alpha,~,~] = quadprog(H,[],[],[],[],[],lb,[],[],options);
abundance=alpha(1:m).*(abs(alpha(1:m))>1e-8);% If less than 1e-8, consider the value as 0
if abs(sum(abundance))>1e-8
    abundance = abundance/sum(abundance);
else
    abundance = zeors(m,1);
end
end