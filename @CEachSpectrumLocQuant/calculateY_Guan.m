function Y=calculateY_Guan(~,vNonRedunTheoryIonMz,matchedExpPeaks)
% Calculate the Y matrix in $Y=X\alpha+\epsilon$
% Input: 
%   vNonRedunTheoryIonMz - Site-discrimining ions, each row is a fragment ion:
%       [m/z, type (1 is b ion, 2 is y ion), ion number (position), charge, 
%       number of modifications, class index, whether an IMP can generate this ion]
%   matchedExpPeaks - experimental spectrum peaks, first column is m/z, second column is normalized intensity, third column is the true intensity obtained experimentally
%   ms2_tolerance - fragment ion matching error
% Output: 
%   Y - the Y matrix in $Y=X\alpha+\epsilon$
Y=zeros(size(vNonRedunTheoryIonMz,1),1);
intenSum=zeros(max(vNonRedunTheoryIonMz(:,6)),1);   % intensity sum for the fragment ion in the same class
for iPeak=1:size(matchedExpPeaks,1)
    intenSum(vNonRedunTheoryIonMz(matchedExpPeaks(iPeak,1),6)) = ...
        intenSum(vNonRedunTheoryIonMz(matchedExpPeaks(iPeak,1),6))+matchedExpPeaks(iPeak,2);
end
Y(matchedExpPeaks(:,1))=matchedExpPeaks(:,2)./intenSum(vNonRedunTheoryIonMz(matchedExpPeaks(:,1),6));
end