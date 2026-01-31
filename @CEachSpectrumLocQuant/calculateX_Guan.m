% Calculate the X matrix in $Y=X\alpha+\epsilon$
% Input: 
%   vNonRedunTheoryIonMz (L x M double) - Site-discrimining ions, each row is a fragment ion:
%       [m/z, type (1 is b ion, 2 is y ion), ion number (position), charge, 
%       number of modifications, class index, whether an IMP can generate this ion]
% Output: 
%   X (L x P double) - the X matrix in $Y=X\alpha+\epsilon$
function [X]=calculateX_Guan(~,vNonRedunTheoryIonMz)
X=vNonRedunTheoryIonMz(:,7:end);
end

