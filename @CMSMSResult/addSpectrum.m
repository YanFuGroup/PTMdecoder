function addSpectrum(obj, datasetName, spectrumName)
    % Input:
    %   datasetName (1 x 1 char/string)
    %       dataset name
    %   spectrumName (1 x 1 char/string)
    %       spectrum name
    if obj.CurrentPeptideIdx == 0
        error('CMSMSResult:NoPeptide', 'Cannot add spectrum without a peptide context.');
    end
    
    obj.CurrentSpectrumIdx = obj.CurrentSpectrumIdx + 1;
    
    % Check capacity and buffer if needed
    % Using direct indexing for speed check
    currentCap = length(obj.Peptides(obj.CurrentPeptideIdx).spectrum_list);
    if obj.CurrentSpectrumIdx > currentCap
        % TODO: buffer size as parameter?
        % Extend buffer by 20
        % Touching the end element automatically expands the struct array with empty fields
        obj.Peptides(obj.CurrentPeptideIdx).spectrum_list(currentCap + 20).peptidoform_num = 0;
    end
    
    % Expand spectrum list for current peptide
    obj.Peptides(obj.CurrentPeptideIdx).spectrum_list(obj.CurrentSpectrumIdx).dataset_name = datasetName;
    obj.Peptides(obj.CurrentPeptideIdx).spectrum_list(obj.CurrentSpectrumIdx).spectrum_name = spectrumName;
    obj.Peptides(obj.CurrentPeptideIdx).spectrum_list(obj.CurrentSpectrumIdx).peptidoform_num = 0;
    
    % Initialize buffers (empty)
    obj.Peptides(obj.CurrentPeptideIdx).spectrum_list(obj.CurrentSpectrumIdx).peptidoform_list_str = {};
    obj.Peptides(obj.CurrentPeptideIdx).spectrum_list(obj.CurrentSpectrumIdx).peptidoform_list_abun = [];
end