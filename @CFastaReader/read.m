function mapProt = read(obj, reject_conflicting_sequences)
% Read sequence from file to mapProt
% Input:
%   obj (CFastaReader)
%       FASTA reader instance
%   reject_conflicting_sequences (logical, optional; default false)
%       Reject different sequences resolving to the same protein name.
%       Site numbering opts in; existing callers retain last-entry-wins.
% Output:
%   mapProt (containers.Map)
%       key: protein name, value: sequence string

if nargin < 2
    reject_conflicting_sequences = false;
end

% Initialize map
mapProt = containers.Map('KeyType', 'char', 'ValueType', 'char');

fin = fopen(obj.m_strFilePath, 'r');
if fin <= 0
    error(['Can not open file: ', obj.m_strFilePath]);
end
cleanup_file = onCleanup(@() fclose(fin)); %#ok<NASGU>

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
            storeSequence(mapProt, currentProtName, currentSeq, reject_conflicting_sequences);
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
    storeSequence(mapProt, currentProtName, currentSeq, reject_conflicting_sequences);
end

end

function storeSequence(mapProt, name, sequence, reject_conflicts)
% Check before assignment so conflicting evidence cannot be overwritten.
if reject_conflicts && isKey(mapProt, name) && ~strcmp(mapProt(name), sequence)
    error('CFastaReader:ConflictingProteinSequence', ...
        'Conflicting FASTA sequences resolve to protein name "%s".', name);
end
mapProt(name) = sequence;
end
