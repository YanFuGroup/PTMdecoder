% Core program for solving the optimization problem, does not depend on individual spectra, constant model of fragmentation efficiency, solved by least squares
% Input: 
%   X - the X matrix in $Y=X\alpha+\epsilon$
%   Y - the Y vector in $Y=X\alpha+\epsilon$
% Output: 
%   abundance - the \alpha vector in $Y=X\alpha+\epsilon$, representing the relative abundance of each IMP
function [abundance]=coreGuan_OLS(~,X,Y)
H=X'*X;
f=-X'*Y;
lb=zeros(size(X,2),1);
options = optimoptions('quadprog','Display','off');
[abundance,~,~] = quadprog(H,f,[],[],[],[],lb,[],[],options);
abundance = abundance.*(abs(abundance)>1e-8);
if abs(sum(abundance))>1e-8
    abundance = abundance/sum(abundance);
else
    abundance = zeros(length(abundance));
end
end