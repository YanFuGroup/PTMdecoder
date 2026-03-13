function addPeptide(obj, sequence)
% Input:
%   sequence (1 x 1 char/string)
%       peptide sequence

obj.CurrentPeptideIdx = length(obj.Peptides) + 1;
obj.CurrentSpectrumIdx = 0;

% Expand structure
obj.Peptides(obj.CurrentPeptideIdx).peptide_sequence = sequence;
obj.Peptides(obj.CurrentPeptideIdx).spectrum_list = ...
    struct('dataset_name',{},'spectrum_name',{}, ...
    'jaccard_stability',{}, ...
    'peptidoform_list_str',{},'peptidoform_list_abun',{}, ...
    'peptidoform_list_support_freq',{},'peptidoform_num',{});
end