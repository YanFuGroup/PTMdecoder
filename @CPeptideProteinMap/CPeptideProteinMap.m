classdef CPeptideProteinMap < handle
    % CPeptideProteinMap establishes a mapping from peptides to proteins.
    % Pass a tab-delimited txt file through the constructor, read the peptide and protein columns.
    
    properties (Access = private)
        m_map % containers.Map storing the mapping, Key is peptide, Value is cell array of proteins
    end
    
    methods
        function obj = CPeptideProteinMap(filePath)
            % CPeptideProteinMap constructor
            % Input:
            %   filePath (1 x 1 char/string)
            %       Path to the tab-delimited text file
            
            % Initialize Map
            obj.m_map = containers.Map('KeyType', 'char', 'ValueType', 'any');
            
            if nargin > 0
                obj.load_file(filePath);
            end
        end
        
        % Get proteins for a peptide
        proteinCell = get_proteins(obj, peptide)
    end
    
    methods (Access = private)
        % Load file and parse data
        load_file(obj, filePath)
    end
end
