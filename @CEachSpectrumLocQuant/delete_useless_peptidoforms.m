function [massArrangement, vNonRedunTheoryIonMz] = delete_useless_peptidoforms( ...
        ~, matchedExpPeaks, massArrangement, vNonRedunTheoryIonMz)
% Delete the peptidoforms that have no matched peaks and the peptidoforms that are subsets of other peptidoforms.
% Input:
%   matchedExpPeaks:
%       List of matched experimental peaks.
%       The first column is the index of the matched vNonRedunTheoryIonMz, the second column is the normalized intensity, the third column is the real intensity obtained from the experiment
%   massArrangement:
%       The mass arrangement of the peptidoforms.
%       Each row represent an IMP, each column corresponds to a modification site, in the order of the site in the peptide sequence
%   vNonRedunTheoryIonMz:
%       Theoretical non-redundant ions.
%       Each row is a fragment ion, each column represents its various basic information, including modifications
%       [m/z of the ion at p charge, type (1 is b ion, 2 is y ion), ion index (position), charge,
%       number of modifications, type index charge group, which IMP can generate this ion]
% Output:
%   massArrangement, vNonRedunTheoryIonMz:

% Get the number of non-redundant ions
num_nrti = size(vNonRedunTheoryIonMz, 2) - 6;
usedNonRedunTheoryIonMz = vNonRedunTheoryIonMz(matchedExpPeaks(:, 1), :);

% Initialize the indices of the peptidoforms that have matched peaks
idx_peptidoforms = [];

% Traverse each peptidoform and delete the peptidoforms that have no matched peaks
for i_form = 1:num_nrti
    if sum(usedNonRedunTheoryIonMz(:, i_form+6)) > 0
        idx_peptidoforms = [idx_peptidoforms, i_form]; %#ok<AGROW> 
    end
end

% Check the relationship between elements in idx_peptidoforms and remove the ones that meet the condition:
%   If the fragment ions produced by one peptidoform are a subset (proper subset) of the fragment ions produced by another peptidoform, then remove the former.
to_delete = false(1, length(idx_peptidoforms));

for i = 1:length(idx_peptidoforms)
    if to_delete(i)
        continue;
    end
    for j = 1:length(idx_peptidoforms)
        if i ~= j && ~to_delete(j)
            diff = usedNonRedunTheoryIonMz(:, idx_peptidoforms(i)+6) - usedNonRedunTheoryIonMz(:, idx_peptidoforms(j)+6);
            if all(diff <= 0) && any(diff < 0)
                to_delete(i) = true; % Mark for deletion
                break;
            end
        end
    end
end

% Remove marked elements from idx_peptidoforms
idx_peptidoforms = idx_peptidoforms(~to_delete);

% Update the mass arrangement and the non-redundant ions
massArrangement = massArrangement(idx_peptidoforms, :);
vNonRedunTheoryIonMz = vNonRedunTheoryIonMz(:, [1:6, idx_peptidoforms+6]);

end