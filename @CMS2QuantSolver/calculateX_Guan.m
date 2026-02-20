function [X]=calculateX_Guan(vNonRedunTheoryIonMz)
% calculateX_Guan - Build X matrix for Guan/equal-FE models.
% Input:
%   vNonRedunTheoryIonMz (L x T double)
%       Non-redundant ion table with columns:
%       [m/z, type (1 is b ion, 2 is y ion), ion number (position), charge, number of modifications, class index, whether an IMP can generate this ion]
% Output:
%   X (L x M double)
%       X matrix for Guan/equal-FE models.
%       Each row is an ion, each column is an IMP.
X=vNonRedunTheoryIonMz(:,7:end);
end