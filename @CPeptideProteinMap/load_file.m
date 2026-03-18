function load_file(obj, filePath)
    % load_file reads and parses the tab-delimited file
    % Inputs:
    %   obj (CPeptideProteinMap)
    %       Peptide->protein mapping instance
    %   filePath (1 x 1 char/string)
    %       Path to the tab-delimited text file
    
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
    lastPct = -1;
    for i = 1:numRows
        if numRows > 0
            currentPct = floor(double(i) / double(numRows) * 100);
            if i == 1 || currentPct > lastPct || i == numRows
                CLogger.progress('peptide_protein_map_load', i, numRows);
                lastPct = currentPct;
            end
        end

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

    if numRows > 0 && lastPct < 100
        CLogger.progress('peptide_protein_map_load', numRows, numRows);
    end
end