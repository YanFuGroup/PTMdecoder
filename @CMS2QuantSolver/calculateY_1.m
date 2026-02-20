function [Y]=calculateY_1(vNonRedunTheoryIonMz,Peaks)
% calculateY_1 - Build Y for equal-FE model directly from matched peak intensities.
% Input:
%   vNonRedunTheoryIonMz (L x T double)
%      Non-redundant ion table with columns:
%      [m/z, type (1 is b ion, 2 is y ion), ion number (position), charge, 
%       number of modifications, class index, whether an IMP can generate this ion]
%   Peaks (K x 3 double)
%       Matched peaks [ion_index, normalized_intensity, raw_intensity].
% Output:
%   Y (L x 1 double)

Y=zeros(size(vNonRedunTheoryIonMz,1),1);
Y(Peaks(:,1))=Peaks(:,2);
end
