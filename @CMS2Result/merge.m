function mergedResult = merge(resultList)
% Merge multiple CMS2Result objects while preserving all hierarchy metadata.
% Input:
%   resultList (1 x N cell or array)
%       CMS2Result objects to merge in order
% Output:
%   mergedResult (CMS2Result)
%       merged result with peptides grouped by sequence and all spectra appended
% Notes:
%   - Spectrum and peptidoform entries are copied through public add APIs so
%     the internal peptide index and buffers stay consistent.
%   - Repeated peptide sequences are intentionally merged into one peptide entry.

if nargin < 1 || isempty(resultList)
    mergedResult = CMS2Result();
    return;
end

if ~iscell(resultList)
    resultList = num2cell(resultList);
end

mergedResult = CMS2Result();
for i_result = 1:numel(resultList)
    srcResult = resultList{i_result};
    if ~isa(srcResult, 'CMS2Result')
        error('CMS2Result:InvalidMergeInput', ...
            'merge expects CMS2Result objects.');
    end

    for i_pep = 1:numel(srcResult.Peptides)
        srcPeptide = srcResult.Peptides(i_pep);
        mergedResult.addOrSelectPeptide(srcPeptide.peptide_sequence);

        for i_spec = 1:numel(srcPeptide.spectrum_list)
            srcSpectrum = srcPeptide.spectrum_list(i_spec);
            spectrumMeta = struct( ...
                'jaccard_stability', getStructFieldOrDefault(srcSpectrum, 'jaccard_stability', NaN), ...
                'vif_all_imp_max', getStructFieldOrDefault(srcSpectrum, 'vif_all_imp_max', NaN), ...
                'vif_reported_imp_max', getStructFieldOrDefault(srcSpectrum, 'vif_reported_imp_max', NaN));
            mergedResult.addSpectrum(srcSpectrum.dataset_name, srcSpectrum.spectrum_name, spectrumMeta);

            peptidoformNum = srcSpectrum.peptidoform_num;
            for i_form = 1:peptidoformNum
                peptidoformMeta = struct( ...
                    'support_frequency', getArrayValueOrDefault(srcSpectrum, 'peptidoform_list_support_freq', i_form, NaN), ...
                    'vif', getArrayValueOrDefault(srcSpectrum, 'peptidoform_list_vif', i_form, NaN), ...
                    'abundance_mad', getArrayValueOrDefault(srcSpectrum, 'peptidoform_list_abundance_mad', i_form, NaN));
                mergedResult.addPeptidoform( ...
                    srcSpectrum.peptidoform_list_str{i_form}, ...
                    srcSpectrum.peptidoform_list_abun(i_form), ...
                    peptidoformMeta);
            end
        end
    end
end

mergedResult.compress();
end


function value = getStructFieldOrDefault(sourceStruct, fieldName, defaultValue)
% Return a struct field when present, otherwise a caller-provided default.
if isfield(sourceStruct, fieldName)
    value = sourceStruct.(fieldName);
else
    value = defaultValue;
end
end


function value = getArrayValueOrDefault(sourceStruct, fieldName, idx, defaultValue)
% Return an indexed metadata value when present and populated.
if isfield(sourceStruct, fieldName) && numel(sourceStruct.(fieldName)) >= idx
    value = sourceStruct.(fieldName)(idx);
else
    value = defaultValue;
end
end
