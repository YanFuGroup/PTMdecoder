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
            obj.Peptides = struct('peptide_sequence', {}, 'spectrum_list', {});
            obj.CurrentPeptideIdx = 0;
        end
        
        function addPeptide(obj, sequence)
            obj.CurrentPeptideIdx = obj.CurrentPeptideIdx + 1;
            obj.CurrentSpectrumIdx = 0;
            
            % Expand structure
            obj.Peptides(obj.CurrentPeptideIdx).peptide_sequence = sequence;
            obj.Peptides(obj.CurrentPeptideIdx).spectrum_list = ...
                struct('dataset_name',{},'spectrum_name',{}, ...
                'peptidoform_list_str',{},'peptidoform_list_abun',{},'peptidoform_num',{});
        end
        
        function addSpectrum(obj, datasetName, spectrumName)
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
        
        function addPeptidoform(obj, peptidoform_str, relative_abundance)
            if obj.CurrentPeptideIdx == 0 || obj.CurrentSpectrumIdx == 0
                error('CMSMSResult:NoSpectrum', 'Cannot add peptidoform without peptide and spectrum context.');
            end
            
            % Quick access references
            currentNum = obj.Peptides(obj.CurrentPeptideIdx).spectrum_list(obj.CurrentSpectrumIdx).peptidoform_num;
            newNum = currentNum + 1;
            
            % Update count
            obj.Peptides(obj.CurrentPeptideIdx).spectrum_list(obj.CurrentSpectrumIdx).peptidoform_num = newNum;
            
            % Buffer check logic from original code
            currentLen = length(obj.Peptides(obj.CurrentPeptideIdx).spectrum_list(obj.CurrentSpectrumIdx).peptidoform_list_abun);
            if newNum > currentLen
                % Extend by 50
                % TODO: make the buffer size a parameter?
                obj.Peptides(obj.CurrentPeptideIdx).spectrum_list(obj.CurrentSpectrumIdx).peptidoform_list_str{newNum + 50} = '';
                obj.Peptides(obj.CurrentPeptideIdx).spectrum_list(obj.CurrentSpectrumIdx).peptidoform_list_abun(newNum + 50) = 0;
            end
            
            obj.Peptides(obj.CurrentPeptideIdx).spectrum_list(obj.CurrentSpectrumIdx).peptidoform_list_str{newNum} = peptidoform_str;
            obj.Peptides(obj.CurrentPeptideIdx).spectrum_list(obj.CurrentSpectrumIdx).peptidoform_list_abun(newNum) = relative_abundance;
        end
        
        function compress(obj)
            % Organize the data structure: delete the empty peptide, spectrum, and peptidoform
            % Also trim the buffers
            
            % Iterate backwards to safely delete
            for i = length(obj.Peptides):-1:1
                if isempty(obj.Peptides(i).spectrum_list)
                    obj.Peptides(i) = [];
                else
                    for j = length(obj.Peptides(i).spectrum_list):-1:1
                        num = obj.Peptides(i).spectrum_list(j).peptidoform_num;
                        % Check for 0 or empty (uninitialized buffer slots)
                        if isempty(num) || num == 0
                            obj.Peptides(i).spectrum_list(j) = [];
                        else
                            % Trim buffers
                             obj.Peptides(i).spectrum_list(j).peptidoform_list_str(num+1:end) = [];
                             obj.Peptides(i).spectrum_list(j).peptidoform_list_abun(num+1:end) = [];
                        end
                    end
                    if isempty(obj.Peptides(i).spectrum_list)
                        obj.Peptides(i) = [];
                    end
                end
            end
            
            % Reset indices just in case
            obj.CurrentPeptideIdx = length(obj.Peptides);
            if obj.CurrentPeptideIdx > 0
                obj.CurrentSpectrumIdx = length(obj.Peptides(end).spectrum_list);
            else
                obj.CurrentSpectrumIdx = 0;
            end
        end
    end
end
