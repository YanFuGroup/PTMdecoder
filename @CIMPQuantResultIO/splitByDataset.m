function output_paths = splitByDataset(input_path, output_dir, output_prefix)
% splitByDataset
% Split a peptide-level IMP quantification result by Dataset/mgf.
% Input:
%   input_path (1 x 1 char/string)
%       combined peptide-level quantification result path
%   output_dir (1 x 1 char/string)
%       directory for split output files
%   output_prefix (1 x 1 char/string)
%       output file prefix, usually the combined output basename
% Output:
%   output_paths (1 x N cell)
%       generated split output file paths

if nargin < 3 || isempty(output_prefix)
    [~, output_prefix, ~] = fileparts(char(string(input_path)));
end

input_path = char(string(input_path));
output_dir = char(string(output_dir));
output_prefix = char(string(output_prefix));

fin = fopen(input_path, 'r');
if fin < 0
    error('CIMPQuantResultIO:CannotOpenSplitInput', ...
        'Cannot open the report file "%s".', input_path);
end
cleanup_input = onCleanup(@() fclose(fin));

header_lines = cell(1, 3);
for idx_header = 1:3
    header_lines{idx_header} = fgetl(fin);
    if ~ischar(header_lines{idx_header})
        error('CIMPQuantResultIO:InvalidSplitInputHeader', ...
            'Expected three header lines in "%s".', input_path);
    end
end

dataset_index_map = containers.Map('KeyType', 'char', 'ValueType', 'double');
dataset_order = {};
dataset_lines = {};
dataset_line_counts = zeros(1, 0);
dataset_line_capacities = zeros(1, 0);
current_protein_line = '';
block_dataset_written = containers.Map('KeyType', 'char', 'ValueType', 'logical');
current_record_dataset = '';
current_record_lines = {};

    function appendCurrentRecord()
        if isempty(current_record_dataset)
            return;
        end

        dataset_name = current_record_dataset;
        if ~dataset_index_map.isKey(dataset_name)
            dataset_order{end + 1} = dataset_name; %#ok<AGROW>
            dataset_index_map(dataset_name) = numel(dataset_order);
            dataset_lines{end + 1} = cell(1, 64); %#ok<AGROW>
            dataset_line_counts(end + 1) = 0; %#ok<AGROW>
            dataset_line_capacities(end + 1) = 64; %#ok<AGROW>
            appendLinesToDataset(numel(dataset_order), header_lines);
        end

        idx_dataset_current = dataset_index_map(dataset_name);
        if ~block_dataset_written.isKey(dataset_name)
            appendLinesToDataset(idx_dataset_current, {current_protein_line});
            block_dataset_written(dataset_name) = true;
        end
        appendLinesToDataset(idx_dataset_current, current_record_lines);

        current_record_dataset = '';
        current_record_lines = {};
    end

    function appendLinesToDataset(idx_dataset_append, new_lines)
        needed_count = dataset_line_counts(idx_dataset_append) + numel(new_lines);
        if needed_count > dataset_line_capacities(idx_dataset_append)
            new_capacity = max(needed_count, dataset_line_capacities(idx_dataset_append) * 2);
            lines = dataset_lines{idx_dataset_append};
            lines(end + 1:new_capacity) = {''};
            dataset_lines{idx_dataset_append} = lines;
            dataset_line_capacities(idx_dataset_append) = new_capacity;
        end

        idx_start = dataset_line_counts(idx_dataset_append) + 1;
        dataset_lines{idx_dataset_append}(idx_start:needed_count) = new_lines;
        dataset_line_counts(idx_dataset_append) = needed_count;
    end

while ~feof(fin)
    line = fgetl(fin);
    if ~ischar(line) || isempty(line)
        continue;
    end

    switch line(1)
        case '*'
            appendCurrentRecord();
            current_record_dataset = parseDatasetName(line);
            current_record_lines = {line};
        case '@'
            if ~isempty(current_record_dataset)
                current_record_lines{end + 1} = line; %#ok<AGROW>
            end
        otherwise
            appendCurrentRecord();
            current_protein_line = line;
            block_dataset_written = containers.Map('KeyType', 'char', 'ValueType', 'logical');
    end
end
appendCurrentRecord();
clear cleanup_input;

CPathResolver.ensureDir(output_dir);
output_paths = cell(1, numel(dataset_order));
safe_name_dataset_map = containers.Map('KeyType', 'char', 'ValueType', 'char');
for idx_dataset = 1:numel(dataset_order)
    dataset_name = dataset_order{idx_dataset};
    safe_dataset_name = makeSafeDatasetName(dataset_name);
    safe_name_dataset_map = assertUniqueDatasetName( ...
        safe_dataset_name, dataset_name, safe_name_dataset_map);
    output_path = fullfile(output_dir, [output_prefix, '__', safe_dataset_name, '.txt']);
    lines = dataset_lines{idx_dataset}(1:dataset_line_counts(idx_dataset));
    writeLines(output_path, lines);
    output_paths{idx_dataset} = output_path;
end
end


function dataset_name = parseDatasetName(record_line)
% Parse the Dataset/mgf column from a peptide quant record line.
parts = regexp(record_line, '\t', 'split');
if numel(parts) < 4 || ~strcmp(parts{1}, '*')
    error('CIMPQuantResultIO:InvalidSplitRecordLine', ...
        'Record line must contain a Dataset column: %s', record_line);
end
dataset_name = parts{4};
end


function safe_name = makeSafeDatasetName(dataset_name)
% Convert a Dataset/mgf value into a stable filesystem-safe basename.
dataset_name = char(string(dataset_name));
dataset_name = strrep(dataset_name, '\', '/');
slash_idx = find(dataset_name == '/', 1, 'last');
if ~isempty(slash_idx)
    dataset_name = dataset_name(slash_idx + 1:end);
end

if length(dataset_name) >= 4 && strcmpi(dataset_name(end - 3:end), '.mgf')
    dataset_name = dataset_name(1:end - 4);
end

safe_name = regexprep(dataset_name, '[^A-Za-z0-9._-]', '_');
if isempty(safe_name)
    safe_name = 'dataset';
end
end


function safe_name_dataset_map = assertUniqueDatasetName(safe_name, dataset_name, safe_name_dataset_map)
% Fail fast when different Dataset values would write the same split file.
if safe_name_dataset_map.isKey(safe_name)
    existing_dataset_name = safe_name_dataset_map(safe_name);
    if ~strcmp(existing_dataset_name, dataset_name)
        error('CIMPQuantResultIO:SplitDatasetNameCollision', ...
            ['Dataset values "%s" and "%s" both map to split output basename "%s". ', ...
            'Rename one Dataset or choose non-conflicting input names.'], ...
            existing_dataset_name, dataset_name, safe_name);
    end
    return;
end

safe_name_dataset_map(safe_name) = dataset_name;
end


function writeLines(output_path, lines)
% Write prepared report lines to one split output file.
fid = fopen(output_path, 'w');
if fid < 0
    error('CIMPQuantResultIO:CannotOpenSplitOutput', ...
        'Cannot open the split report file "%s".', output_path);
end
cleanup_output = onCleanup(@() fclose(fid));

for idx_line = 1:numel(lines)
    fprintf(fid, '%s\n', lines{idx_line});
end
clear cleanup_output;
end
