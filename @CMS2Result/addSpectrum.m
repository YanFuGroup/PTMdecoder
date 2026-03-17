function addSpectrum(obj, datasetName, spectrumName, varargin)
% Input:
%   datasetName (1 x 1 char/string)
%       dataset name
%   spectrumName (1 x 1 char/string)
%       spectrum name
%   varargin{1} (optional, 1 x 1 struct)
%       optional spectrum metadata
%       - jaccard_stability (1 x 1 double)

if obj.CurrentPeptideIdx == 0
    error('CMS2Result:NoPeptide', 'Cannot add spectrum without a peptide context.');
end

jaccard_stability = NaN;
if nargin >= 4 && ~isempty(varargin{1})
    spectrumMeta = varargin{1};
    if ~isstruct(spectrumMeta)
        error('CMS2Result:InvalidSpectrumMetadata', ...
            'Optional spectrum metadata must be a struct.');
    end
    if isfield(spectrumMeta, 'jaccard_stability')
        jaccard_stability = spectrumMeta.jaccard_stability;
    end
end

obj.CurrentSpectrumIdx = obj.CurrentSpectrumIdx + 1;

% Expand spectrum list for current peptide
obj.Peptides(obj.CurrentPeptideIdx).spectrum_list(obj.CurrentSpectrumIdx).dataset_name = datasetName;
obj.Peptides(obj.CurrentPeptideIdx).spectrum_list(obj.CurrentSpectrumIdx).spectrum_name = spectrumName;
obj.Peptides(obj.CurrentPeptideIdx).spectrum_list(obj.CurrentSpectrumIdx).jaccard_stability = jaccard_stability;
obj.Peptides(obj.CurrentPeptideIdx).spectrum_list(obj.CurrentSpectrumIdx).peptidoform_num = 0;

% Initialize buffers (empty)
obj.Peptides(obj.CurrentPeptideIdx).spectrum_list(obj.CurrentSpectrumIdx).peptidoform_list_str = {};
obj.Peptides(obj.CurrentPeptideIdx).spectrum_list(obj.CurrentSpectrumIdx).peptidoform_list_abun = [];
obj.Peptides(obj.CurrentPeptideIdx).spectrum_list(obj.CurrentSpectrumIdx).peptidoform_list_support_freq = [];
obj.Peptides(obj.CurrentPeptideIdx).spectrum_list(obj.CurrentSpectrumIdx).peptidoform_list_abundance_mad = [];
end