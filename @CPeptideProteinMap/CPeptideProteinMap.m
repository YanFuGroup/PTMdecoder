classdef CPeptideProteinMap < handle
    % CPeptideProteinMap establishes a mapping from peptides to proteins.
    % Pass a tab-delimited txt file through the constructor, read the peptide and protein columns.
    
    properties (Access = private)
        m_map % containers.Map storing the mapping, Key is peptide, Value is cell array of proteins
    end
    
    methods
        function obj = CPeptideProteinMap(filePath)
            % CPeptideProteinMap constructor
            % filePath: path to the tab-delimited text file
            
            % Initialize Map
            obj.m_map = containers.Map('KeyType', 'char', 'ValueType', 'any');
            
            if nargin > 0
                obj.load_file(filePath);
            end
        end
        
        function proteinCell = get_proteins(obj, peptide)
            % get_proteins takes a peptide string as input and returns a cell array of corresponding protein names
            % peptide: peptide string
            % proteinCell: cell array of protein names, returns {} if not found
            
            if ischar(peptide)
                pepKey = peptide;
            elseif isstring(peptide)
                pepKey = char(peptide);
            else
                proteinCell = {};
                return;
            end
            
            if isKey(obj.m_map, pepKey)
                proteinCell = obj.m_map(pepKey);
            else
                proteinCell = {};
            end
        end
    end
    
    methods (Access = private)
        function load_file(obj, filePath)
            % load_file reads and parses the tab-delimited file
            
            if ~exist(filePath, 'file')
                error('File not found: %s', filePath);
            end
            
            % Detect import options
            try
                opts = detectImportOptions(filePath, 'FileType', 'text', 'Delimiter', '\t');
                opts.VariableNamingRule = 'preserve';
                
                % Ensure all columns are read as strings to avoid automatic conversion to numbers or categories
                for i = 1:length(opts.VariableNames)
                    opts = setvaropts(opts, opts.VariableNames{i}, 'Type', 'char');
                end
                
                tbl = readtable(filePath, opts);
            catch
                % If detectImportOptions fails, try simple reading
                tbl = readtable(filePath, 'FileType', 'text', 'Delimiter', '\t', 'VariableNamingRule', 'preserve');
            end
            
            % Find column names (case-insensitive)
            colNames = tbl.Properties.VariableNames;
            pepColIdx = find(strcmpi(colNames, 'peptide'), 1);
            protColIdx = find(strcmpi(colNames, 'protein'), 1);
            
            if isempty(pepColIdx) || isempty(protColIdx)
                error('Required ''peptide'' or ''protein'' columns not found in the file.');
            end
            
            % Iterate through each row to populate the mapping
            numRows = height(tbl);
            for i = 1:numRows
                % Get the content from the original data cells
                pepCell = tbl{i, pepColIdx};
                protCell = tbl{i, protColIdx};
                
                if iscell(pepCell)
                    pepStr = char(pepCell{1});
                else
                    pepStr = char(pepCell);
                end
                
                if iscell(protCell)
                    protStr = char(protCell{1});
                else
                    protStr = char(protCell);
                end
                
                % Trim leading and trailing whitespace
                pepStr = strtrim(pepStr);
                protStr = strtrim(protStr);
                
                if isempty(pepStr)
                    continue;
                end
                
                % Parse proteins (comma-separated)
                protList = strsplit(protStr, ',');
                protList = strtrim(protList);
                protList = protList(~cellfun(@isempty, protList));
                
                % Store in Map, supporting merging for peptides appearing in multiple rows
                if isKey(obj.m_map, pepStr)
                    existingProts = obj.m_map(pepStr);
                    obj.m_map(pepStr) = unique([existingProts, protList], 'stable');
                else
                    obj.m_map(pepStr) = protList;
                end
            end
        end
    end
end
