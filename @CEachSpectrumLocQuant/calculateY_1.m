function [Y]=calculateY_1(~,vNonRedunTheoryIonMz,Peaks)
% Calculate the Y matrix in $Y=X\alpha+\epsilon$
% Input: 
%   vNonRedunTheoryIonMz - Site-discrimining ions, each row is a fragment ion:
%       [m/z, type (1 is b ion, 2 is y ion), ion number (position), charge, 
%       number of modifications, class index, whether an IMP can generate this ion]
%   Peaks - experimental spectrum peaks, first column is m/z, second column is normalized intensity, third column is the true intensity obtained experimentally
%   ms2_tolerance - fragment ion matching error
% Output: 
%   Y - the Y matrix in $Y=X\alpha+\epsilon$
Y=zeros(size(vNonRedunTheoryIonMz,1),1);
Y(Peaks(:,1))=Peaks(:,2);
end
