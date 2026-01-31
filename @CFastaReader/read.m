function mapProt = read(obj)
% Read sequence from file to mapProt
% Input:
%   obj (CFastaReader)
%       FASTA reader instance
% Output:
%   mapProt (containers.Map)
%       key: protein name, value: sequence string

% Initialize map
mapProt = containers.Map('KeyType', 'char', 'ValueType', 'char');

fin = fopen(obj.m_strFilePath, 'r');
if fin <= 0
    error(['Can not open file: ', obj.m_strFilePath]);
end

currentProtName = '';
currentSeq = '';

while ~feof(fin)
    strLine = strtrim(fgetl(fin));
    if isempty(strLine)
        continue;
    end
    
    if strLine(1) == '>'
        % Save previous protein if it exists
        if ~isempty(currentProtName)
            mapProt(currentProtName) = currentSeq;
        end
        
        % Extract protein name using regular expression
        match = regexp(strLine, obj.m_regular_exp, 'tokens', 'once');
        if ~isempty(match)
            currentProtName = match{1};
        else
            warning(['Protein name not found with regular expression in header: ', strLine]);
            currentProtName = strLine(2:end); % Fallback
        end
        currentSeq = '';
    else
        % Concatenate sequence lines, removing any whitespace
        currentSeq = [currentSeq, regexprep(strLine, '\s+', '')]; %#ok<AGROW>
    end
end

% Save the last protein
if ~isempty(currentProtName)
    mapProt(currentProtName) = currentSeq;
end

fclose(fin);
end
