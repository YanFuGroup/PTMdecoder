function addSpectrum(obj, datasetName, spectrumName)
% Input:
%   datasetName (1 x 1 char/string)
%       dataset name
%   spectrumName (1 x 1 char/string)
%       spectrum name

if obj.CurrentPeptideIdx == 0
    error('CMS2Result:NoPeptide', 'Cannot add spectrum without a peptide context.');
end

obj.CurrentSpectrumIdx = obj.CurrentSpectrumIdx + 1;

% Expand spectrum list for current peptide
obj.Peptides(obj.CurrentPeptideIdx).spectrum_list(obj.CurrentSpectrumIdx).dataset_name = datasetName;
obj.Peptides(obj.CurrentPeptideIdx).spectrum_list(obj.CurrentSpectrumIdx).spectrum_name = spectrumName;
obj.Peptides(obj.CurrentPeptideIdx).spectrum_list(obj.CurrentSpectrumIdx).peptidoform_num = 0;

% Initialize buffers (empty)
obj.Peptides(obj.CurrentPeptideIdx).spectrum_list(obj.CurrentSpectrumIdx).peptidoform_list_str = {};
obj.Peptides(obj.CurrentPeptideIdx).spectrum_list(obj.CurrentSpectrumIdx).peptidoform_list_abun = [];
end