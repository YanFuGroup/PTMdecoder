classdef CMS2Result < handle
    % CMS2Result: Data Object for MS/MS level results
    % Hierarchical structure: Peptide -> Spectrum -> Peptidoform

    properties
        % The hierarchical structure of the overall MS/MS results:
        % struct array with fields:
        %   - peptide_sequence: char/string
        %   - spectrum_list: struct array
        %       - dataset_name: char/string
        %       - spectrum_name: char/string
        %       - jaccard_stability: 1 x 1 double
        %       - peptidoform_list_str: cell array of strings (N x 1)
        %       - peptidoform_list_abun: double array (N x 1)
        %       - peptidoform_list_support_freq: double array (N x 1)
        %       - peptidoform_list_abundance_mad: double array (N x 1)
        %       - peptidoform_num: int
        Peptides
    end

    properties(Access = private)
        CurrentPeptideIdx = 0;
        % TODO: Current spectrum index within the current peptide, not global
        CurrentSpectrumIdx = 0;
    end

    methods
        function obj = CMS2Result()
            % Constructor
            obj.Peptides = struct('peptide_sequence', {}, 'spectrum_list', {});
            obj.CurrentPeptideIdx = 0;
        end

        addPeptide(obj, sequence)
        addOrSelectPeptide(obj, sequence)
        addSpectrum(obj, datasetName, spectrumName, varargin)
        addPeptidoform(obj, peptidoform_str, relative_abundance, varargin)
        compress(obj)
    end
end
