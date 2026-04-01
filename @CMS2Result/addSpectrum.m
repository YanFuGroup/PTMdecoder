function addSpectrum(obj, datasetName, spectrumName, varargin)
% Input:
%   datasetName (1 x 1 char/string)
%       dataset name
%   spectrumName (1 x 1 char/string)
%       spectrum name
%   varargin{1} (optional, 1 x 1 struct)
%       optional spectrum metadata
%       - jaccard_stability (1 x 1 double)
%       - vif_all_imp_max (1 x 1 double)
%       - vif_reported_imp_max (1 x 1 double)

if obj.CurrentPeptideIdx == 0
    error('CMS2Result:NoPeptide', 'Cannot add spectrum without a peptide context.');
end

jaccard_stability = NaN;
vif_all_imp_max = NaN;
vif_reported_imp_max = NaN;
if nargin >= 4 && ~isempty(varargin{1})
    spectrumMeta = varargin{1};
    if ~isstruct(spectrumMeta)
        error('CMS2Result:InvalidSpectrumMetadata', ...
            'Optional spectrum metadata must be a struct.');
    end
    if isfield(spectrumMeta, 'jaccard_stability')
        jaccard_stability = spectrumMeta.jaccard_stability;
    end
    if isfield(spectrumMeta, 'vif_all_imp_max')
        vif_all_imp_max = spectrumMeta.vif_all_imp_max;
    end
    if isfield(spectrumMeta, 'vif_reported_imp_max')
        vif_reported_imp_max = spectrumMeta.vif_reported_imp_max;
    end
end

currentSpectrumIdx = length(obj.Peptides(obj.CurrentPeptideIdx).spectrum_list) + 1;

% Expand spectrum list for current peptide
obj.Peptides(obj.CurrentPeptideIdx).spectrum_list(currentSpectrumIdx).dataset_name = datasetName;
obj.Peptides(obj.CurrentPeptideIdx).spectrum_list(currentSpectrumIdx).spectrum_name = spectrumName;
obj.Peptides(obj.CurrentPeptideIdx).spectrum_list(currentSpectrumIdx).jaccard_stability = jaccard_stability;
obj.Peptides(obj.CurrentPeptideIdx).spectrum_list(currentSpectrumIdx).vif_all_imp_max = vif_all_imp_max;
obj.Peptides(obj.CurrentPeptideIdx).spectrum_list(currentSpectrumIdx).vif_reported_imp_max = vif_reported_imp_max;
obj.Peptides(obj.CurrentPeptideIdx).spectrum_list(currentSpectrumIdx).peptidoform_num = 0;

% Initialize buffers (empty)
obj.Peptides(obj.CurrentPeptideIdx).spectrum_list(currentSpectrumIdx).peptidoform_list_str = {};
obj.Peptides(obj.CurrentPeptideIdx).spectrum_list(currentSpectrumIdx).peptidoform_list_abun = [];
obj.Peptides(obj.CurrentPeptideIdx).spectrum_list(currentSpectrumIdx).peptidoform_list_support_freq = [];
obj.Peptides(obj.CurrentPeptideIdx).spectrum_list(currentSpectrumIdx).peptidoform_list_vif = [];
obj.Peptides(obj.CurrentPeptideIdx).spectrum_list(currentSpectrumIdx).peptidoform_list_abundance_mad = [];
end