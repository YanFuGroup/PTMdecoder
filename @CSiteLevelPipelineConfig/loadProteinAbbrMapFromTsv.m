function protein_name_abbr_map = loadProteinAbbrMapFromTsv(file_path, col_protein_name, col_abbr_name)
% Load protein-to-abbreviation mapping from a TSV file
% Inputs:
%    file_path (1 x 1 char/string)
%        Path of the TSV mapping file.
%    col_protein_name (1 x 1 char/string)
%        Header column name for protein names.
%    col_abbr_name (1 x 1 char/string)
%        Header column name for abbreviation names.
% Outputs:
%    protein_name_abbr_map (containers.Map)
%        Mapping from normalized protein name to abbreviation.
file_path = normalizeTextInput(file_path, 'file_path');
col_protein_name = strtrim(normalizeTextInput(col_protein_name, 'col_protein_name'));
col_abbr_name = strtrim(normalizeTextInput(col_abbr_name, 'col_abbr_name'));

if isempty(file_path)
    error('CSiteLevelPipelineConfig:InvalidProteinAbbrMapFilePath', ...
        'file_path must not be empty.');
end
if isempty(col_protein_name)
    error('CSiteLevelPipelineConfig:InvalidProteinNameColumn', ...
        'col_protein_name must not be empty.');
end
if isempty(col_abbr_name)
    error('CSiteLevelPipelineConfig:InvalidAbbrColumn', ...
        'col_abbr_name must not be empty.');
end
if ~exist(file_path, 'file')
    error('CSiteLevelPipelineConfig:ProteinAbbrMapFileNotFound', ...
        'Protein abbreviation map file does not exist: %s', file_path);
end

protein_name_abbr_map = containers.Map('KeyType', 'char', 'ValueType', 'char');

fid = fopen(file_path, 'r');
if fid < 0
    error('CSiteLevelPipelineConfig:ProteinAbbrMapFileOpenFailed', ...
        'Cannot open protein abbreviation map file: %s', file_path);
end
cleanup_obj = onCleanup(@() fclose(fid));

header_line = fgetl(fid);
if ~ischar(header_line)
    error('CSiteLevelPipelineConfig:ProteinAbbrMapEmptyFile', ...
        'Protein abbreviation map file is empty: %s', file_path);
end

header_tokens = strsplit(header_line, '	', 'CollapseDelimiters', false);
if isempty(header_tokens)
    error('CSiteLevelPipelineConfig:ProteinAbbrMapInvalidHeader', ...
        'Protein abbreviation map file has invalid header: %s', file_path);
end

header_tokens{1} = stripBom(header_tokens{1});
header_tokens = cellfun(@strtrim, header_tokens, 'UniformOutput', false);

idx_protein = find(strcmp(header_tokens, col_protein_name), 1, 'first');
idx_abbr = find(strcmp(header_tokens, col_abbr_name), 1, 'first');
if isempty(idx_protein)
    error('CSiteLevelPipelineConfig:ProteinAbbrMapMissingProteinColumn', ...
        'Column "%s" is not found in header of file: %s', col_protein_name, file_path);
end
if isempty(idx_abbr)
    error('CSiteLevelPipelineConfig:ProteinAbbrMapMissingAbbrColumn', ...
        'Column "%s" is not found in header of file: %s', col_abbr_name, file_path);
end

num_columns = numel(header_tokens);
format_spec = repmat('%s', 1, num_columns);
column_cells = textscan(fid, format_spec, ...
    'Delimiter', '\t', ...
    'EndOfLine', '\n', ...
    'Whitespace', '', ...
    'MultipleDelimsAsOne', false, ...
    'ReturnOnError', false);

if isempty(column_cells) || isempty(column_cells{1})
    data_row_count = 0;
else
    data_row_count = numel(column_cells{1});
end

protein_col = column_cells{idx_protein};
abbr_col = column_cells{idx_abbr};

protein_col = cellfun(@strtrim, protein_col, 'UniformOutput', false);
abbr_col = cellfun(@strtrim, abbr_col, 'UniformOutput', false);

protein_tokens_by_row = cellfun(@splitProteinTokens, protein_col, 'UniformOutput', false);
abbr_token_by_row = cellfun(@extractFirstNonEmptyAbbrToken, abbr_col, 'UniformOutput', false);
is_valid_row = ~cellfun(@isempty, protein_tokens_by_row) & ~cellfun(@isempty, abbr_token_by_row);
invalid_row_count = 0;

for idx_row = 1:data_row_count
    line_no = idx_row + 1;
    if ~is_valid_row(idx_row)
        invalid_row_count = invalid_row_count + 1;
        continue;
    end

    protein_tokens = protein_tokens_by_row{idx_row};
    abbr_token = abbr_token_by_row{idx_row};

    for idx = 1:numel(protein_tokens)
        protein_key = protein_tokens{idx};

        if protein_name_abbr_map.isKey(protein_key)
            existing_abbr = protein_name_abbr_map(protein_key);
            if ~strcmp(existing_abbr, abbr_token)
                CLogger.error(['[CSiteLevelPipelineConfig:ProteinAbbrMapConflict] ', ...
                    'Conflict at line %d for protein "%s": old="%s", new="%s".'], ...
                    line_no, protein_key, existing_abbr, abbr_token);
            end
            continue;
        end

        protein_name_abbr_map(protein_key) = abbr_token;
    end
end

CLogger.info(['[CSiteLevelPipelineConfig:loadProteinAbbrMapFromTsv] ', ...
    'Loaded %d protein-abbreviation mappings from %s. invalid_rows=%d.'], ...
    protein_name_abbr_map.Count, file_path, invalid_row_count);
end


function tokens = splitProteinTokens(raw_value)
% Split protein field by ';' and keep non-empty trimmed tokens
% Inputs:
%    raw_value (1 x N char/string)
%        Raw protein-name field.
% Outputs:
%    tokens (1 x K cell)
%        Non-empty protein tokens.
parts = strsplit(raw_value, ';');
parts = cellfun(@strtrim, parts, 'UniformOutput', false);
tokens = parts(~cellfun(@isempty, parts));
end


function token = extractFirstNonEmptyAbbrToken(raw_value)
% Extract the first non-empty abbreviation token from ';'-delimited field
% Inputs:
%    raw_value (1 x N char/string)
%        Raw abbreviation field.
% Outputs:
%    token (1 x N char)
%        First non-empty abbreviation token, or empty char.
parts = strsplit(raw_value, ';');
parts = cellfun(@strtrim, parts, 'UniformOutput', false);
first_non_empty_idx = find(~cellfun(@isempty, parts), 1, 'first');
if isempty(first_non_empty_idx)
    token = '';
else
    token = parts{first_non_empty_idx};
end
end


function out = stripBom(in)
% Remove UTF-8 BOM marker from string head
% Inputs:
%    in (1 x N char/string)
%        Input text.
% Outputs:
%    out (1 x N char)
%        Text after BOM removal.
out = char(in);
if ~isempty(out) && double(out(1)) == 65279
    out = out(2:end);
end
end


function out = normalizeTextInput(in, var_name)
% Normalize string-like input into char scalar
% Inputs:
%    in (any)
%        Input value.
%    var_name (1 x N char/string)
%        Variable name for error reporting.
% Outputs:
%    out (1 x N char)
%        Converted char value.
if nargin < 2
    var_name = 'input';
end
if ischar(in)
    out = in;
    return;
end
if isstring(in) && isscalar(in)
    out = char(in);
    return;
end
if isempty(in)
    out = '';
    return;
end
CLogger.error('[CSiteLevelPipelineConfig:InvalidTextType]', ...
    'Expected %s to be char/string scalar.', char(var_name));
end
