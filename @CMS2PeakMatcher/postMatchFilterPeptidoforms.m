function [massArrangement, vNonRedunTheoryIonMz] = postMatchFilterPeptidoforms(matchedExpPeaks, massArrangement, vNonRedunTheoryIonMz)
% postMatchFilterPeptidoforms - Remove peptidoforms with no matched peaks or strict subset ions
% Inputs:
%   matchedExpPeaks (K x 3 double):
%       List of matched experimental peaks.
%       The first column is the index of the matched vNonRedunTheoryIonMz, 
%       the second column is the normalized intensity, 
%       the third column is the real intensity obtained from the experiment
%   massArrangement (M x S double):
%       The mass arrangement of the peptidoforms.
%       Each row represent an IMP, each column corresponds to a modification site, in the order of the site in the peptide sequence
%   vNonRedunTheoryIonMz (L x T double):
%       Theoretical non-redundant ions.
%       Each row is a fragment ion, each column represents its various basic information, including modifications
%       [m/z of the ion at p charge, type (1 is b ion, 2 is y ion), ion index (position), charge,
%       number of modifications, type index charge group, which IMP can generate this ion]
% Outputs:
%   massArrangement (M' x S double):
%      Filtered mass arrangement of the peptidoforms, only including those with matched peaks and not strict subsets of others.
%   vNonRedunTheoryIonMz (L x T' double):
%      Filtered theoretical non-redundant ions, only including those corresponding to the remaining peptidoforms.

% Get the number of non-redundant ions
num_nrti = size(vNonRedunTheoryIonMz, 2) - 6;
usedNonRedunTheoryIonMz = vNonRedunTheoryIonMz(matchedExpPeaks(:, 1), :);

% Initialize the indices of the peptidoforms that have matched peaks
col_sums = sum(usedNonRedunTheoryIonMz(:, 7:6+num_nrti), 1);
IMP_idxs = find(col_sums > 0);

% Check the relationship between elements in IMP_idxs and remove the ones that meet the condition:
%   If the fragment ions produced by one peptidoform are a subset (proper subset) of the fragment ions produced by another peptidoform, then remove the former.
to_delete = false(1, length(IMP_idxs));
for i = 1:length(IMP_idxs)
    if to_delete(i)
        continue;
    end
    for j = 1:length(IMP_idxs)
        if i ~= j && ~to_delete(j)
            diff_val = usedNonRedunTheoryIonMz(:, IMP_idxs(i)+6) - usedNonRedunTheoryIonMz(:, IMP_idxs(j)+6);
            if all(diff_val <= 0) && any(diff_val < 0)
                to_delete(i) = true;
                break;
            end
        end
    end
end

IMP_idxs = IMP_idxs(~to_delete);
massArrangement = massArrangement(IMP_idxs, :);
vNonRedunTheoryIonMz = vNonRedunTheoryIonMz(:, [1:6, IMP_idxs+6]);
end
