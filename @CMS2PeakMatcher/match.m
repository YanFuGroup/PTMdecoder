function matchedExpPeaks = match(expPeaks, vNonRedunTheoryIonMz, ms2_tolerance)
% Iterate through various modifications to find the union of experimental spectrum peak sets that can match the mass-to-charge ratio
% Input: 
%   expPeaks (N x 2 double)
%       experimental spectrum peaks [m/z, intensity]
%   vNonRedunTheoryIonMz (L x 1 double or L x M double)
%       ion info matrix, first column is m/z
% Output: 
%   vMatchedExpPeaks (K x 2 double)
%       matched [ion_index, intensity]
%       ion_index is the row index in vNonRedunTheoryIonMz
% match - Direct m/z matching from experimental peaks to theoretical ions

matchedExpPeaks = [(1:size(vNonRedunTheoryIonMz,1))', zeros(size(vNonRedunTheoryIonMz,1),1)];
for idxExpPeak = 1:size(expPeaks,1)
    iMatchedNonRedun = find(abs(vNonRedunTheoryIonMz(:,1)-expPeaks(idxExpPeak,1)) <= ms2_tolerance);
    for idxMatched = 1:length(iMatchedNonRedun)
        matchedExpPeaks(iMatchedNonRedun(idxMatched),2) = matchedExpPeaks(iMatchedNonRedun(idxMatched),2) + expPeaks(idxExpPeak,2);
    end
end
matchedExpPeaks(matchedExpPeaks(:,2)==0,:) = [];
end
