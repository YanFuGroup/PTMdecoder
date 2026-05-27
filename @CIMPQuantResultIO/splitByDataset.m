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

dataset_lines_map = containers.Map('KeyType', 'char', 'ValueType', 'any');
dataset_order = {};
current_protein_line = '';
block_dataset_written = containers.Map('KeyType', 'char', 'ValueType', 'logical');
current_record_dataset = '';
current_record_lines = {};

    function appendCurrentRecord()
        if isempty(current_record_dataset)
            return;
        end

        dataset_name = current_record_dataset;
        if ~dataset_lines_map.isKey(dataset_name)
            dataset_lines_map(dataset_name) = header_lines;
            dataset_order{end + 1} = dataset_name; %#ok<AGROW>
        end

        lines = dataset_lines_map(dataset_name);
        if ~block_dataset_written.isKey(dataset_name)
            lines{end + 1} = current_protein_line; %#ok<AGROW>
            block_dataset_written(dataset_name) = true;
        end
        lines = [lines, current_record_lines]; %#ok<AGROW>
        dataset_lines_map(dataset_name) = lines;

        current_record_dataset = '';
        current_record_lines = {};
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
for idx_dataset = 1:numel(dataset_order)
    dataset_name = dataset_order{idx_dataset};
    safe_dataset_name = makeSafeDatasetName(dataset_name);
    output_path = fullfile(output_dir, [output_prefix, '__', safe_dataset_name, '.txt']);
    writeLines(output_path, dataset_lines_map(dataset_name));
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
