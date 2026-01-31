classdef CPepNormalization < handle
% Calculate the normalization peptide for sample comparison
% This class encapsulates the functionality from Demo_MCF7_JIB04_normalized_pep.m
% Author:
%   Zhiyuan Cheng
% Create date:
%   2025/07/06
    
    properties (Access = private)
        peptide_list       % List of peptide sequences
        prot_list          % List of corresponding protein names
        result_dir         % Directory for storing results
        spectra_dir        % Directory containing spectral data
        experimentNames    % List of experiment names
        ms1_tolerance      % MS1 mass accuracy tolerance
        resFilterThres     % Threshold for filtering results by relative intensity
        alpha              % Threshold for filtering noise peaks
        input_file_name    % Name of the input file (optional)
        output_file_name   % Name of the output file (optional)
    end
    
    methods
        function obj = CPepNormalization(peptide_list, prot_list, result_dir, ...
                spectra_dir, experimentNames, ms1_tolerance, resFilterThres, alpha, ...
                input_file_name, output_file_name)
            % Constructor
            % Input:
            %   peptide_list (K x 1 cellstr/string)
            %       List of peptide sequences
            %   prot_list (K x 1 cellstr/string)
            %       List of corresponding protein names
            %   result_dir (1 x 1 char/string)
            %       Directory for storing results
            %   spectra_dir (1 x 1 char/string)
            %       Directory containing spectral data
            %   experimentNames (E x 1 cellstr/string)
            %       List of experiment names
            %   ms1_tolerance (struct)
            %       MS1 mass accuracy tolerance structure (fields: value, isppm)
            %   resFilterThres (1 x 1 double)
            %       Threshold for filtering results
            %   alpha (1 x 1 double)
            %       Threshold for filtering noise peaks
            %   input_file_name (1 x 1 char/string, optional)
            %       Name of the input file
            %   output_file_name (1 x 1 char/string, optional)
            %       Name of the output file
            
            obj.peptide_list = peptide_list;
            obj.prot_list = prot_list;
            obj.result_dir = result_dir;
            obj.spectra_dir = spectra_dir;
            obj.experimentNames = experimentNames;
            obj.ms1_tolerance = ms1_tolerance;
            obj.resFilterThres = resFilterThres;
            obj.alpha = alpha;
            if nargin >= 9
                obj.input_file_name = input_file_name;
            else
                obj.input_file_name = 'filtered_result_mascot.txt';
            end
            if nargin >= 10
                obj.output_file_name = output_file_name;
            else
                obj.output_file_name = 'peptide4normalization.txt';
            end

            % Validate input parameters
            obj.validateInputs();
        end
        
        % Run the normalization peptide calculation for all experiments
        runQuant(obj);
    end
    
    methods (Access = private)
        function validateInputs(obj)
            % Validate input parameters
            
            if length(obj.peptide_list) ~= length(obj.prot_list)
                error('peptide_list and prot_list must have the same length');
            end
            
            if ~exist(obj.result_dir, 'dir')
                error('Result directory does not exist: %s', obj.result_dir);
            end
            
            if ~exist(obj.spectra_dir, 'dir')
                error('Spectra directory does not exist: %s', obj.spectra_dir);
            end
            
            if ~isfield(obj.ms1_tolerance, 'value') || ~isfield(obj.ms1_tolerance, 'isppm')
                error('ms1_tolerance must be a struct with fields: value and isppm');
            end
        end
        
        % Process one experiment for normalization peptide calculation
        processOneExperiment(obj, experimentName)
        
        % Create and initialize the output file
        createOutputFile(obj, outputPath);
        
        % Initialize quantification objects for each peptide
        pep_quant = initializeQuantificationObjects(obj, outputPath, ms12DatasetIO);
        
        % Process the filtered results file and extract peptide information
        pep_quant = readSearchResult(obj, fin, input_file_path, ms12DatasetIO, pep_quant);
        
        % Calculate the mass of a peptide sequence
        lfMass = getPeptideMass(obj, pep_seq);
    end
end