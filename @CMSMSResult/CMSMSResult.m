classdef CMSMSResult < handle
    % CMSMSResult: Data Object for MSMS level results
    % Hierarchical structure: Peptide -> Spectrum -> Peptidoform
    
    properties
        % The hierarchical structure of the overall MSMS results: 
        % struct array with fields:
        %   - peptide_sequence: char/string
        %   - spectrum_list: struct array
        %       - dataset_name: char/string
        %       - spectrum_name: char/string
        %       - peptidoform_list_str: cell array of strings (N x 1)
        %       - peptidoform_list_abun: double array (N x 1)
        %       - peptidoform_num: int
        Peptides
    end
    
    properties(Access = private)
        CurrentPeptideIdx = 0;
        % TODO: Current spectrum index within the current peptide, not global
        CurrentSpectrumIdx = 0;
    end
    
    methods
        function obj = CMSMSResult()
            % Constructor
            obj.Peptides = struct('peptide_sequence', {}, 'spectrum_list', {});
            obj.CurrentPeptideIdx = 0;
        end
        
        addPeptide(obj, sequence)
        addSpectrum(obj, datasetName, spectrumName)
        addPeptidoform(obj, peptidoform_str, relative_abundance)
        compress(obj)
    end
end
