function Y=calculateY_Guan(vNonRedunTheoryIonMz,matchedExpPeaks)
% calculateY_Guan - Build Y for Guan model using within-ion-group normalization.
% Input:
%   vNonRedunTheoryIonMz (L x T double)
%      Non-redundant ion table with columns:
%      [m/z, type (1 is b ion, 2 is y ion), ion number (position), charge, 
%       number of modifications, class index, whether an IMP can generate this ion]
%   matchedExpPeaks (K x 3 double)
%       Experimental spectrum peaks, each row is a matched peak:
%       [ion_index, normalized_intensity, raw_intensity].
% Output:
%   Y (L x 1 double)

Y=zeros(size(vNonRedunTheoryIonMz,1),1);
intenSum=zeros(max(vNonRedunTheoryIonMz(:,6)),1);
for iPeak=1:size(matchedExpPeaks,1)
	intenSum(vNonRedunTheoryIonMz(matchedExpPeaks(iPeak,1),6)) = ...
		intenSum(vNonRedunTheoryIonMz(matchedExpPeaks(iPeak,1),6))+matchedExpPeaks(iPeak,2);
end
Y(matchedExpPeaks(:,1))=matchedExpPeaks(:,2)./intenSum(vNonRedunTheoryIonMz(matchedExpPeaks(:,1),6));
end