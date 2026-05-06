function dataset_charge_map = get_dataset_charge_map(obj, dataset_name)
%GET_DATASET_CHARGE_MAP Read all spectrum TITLE/CHARGE pairs from one MGF file.
% Input:
%   obj (CMgfDatasetIO)
%       dataset IO instance
%   dataset_name (1 x 1 char/string)
%       MGF dataset filename
% Output:
%   dataset_charge_map (containers.Map)
%       map from spectrum TITLE or parsed scan number to precursor charge

dataset_name = char(dataset_name);
cache_path = get_charge_cache_path(obj.m_strFoldname, dataset_name);
mgf_info = dir(fullfile(obj.m_strFoldname, dataset_name));
if isempty(mgf_info)
    error('Cannot find MGF dataset: %s', fullfile(obj.m_strFoldname, dataset_name));
end

if isfile(cache_path)
    cache_data = load(cache_path);
    if isfield(cache_data, 'dataset_charge_map') && isfield(cache_data, 'cache_meta') ...
            && is_charge_cache_valid(cache_data.cache_meta, mgf_info)
        dataset_charge_map = cache_data.dataset_charge_map;
        CLogger.debug('Loaded cached MGF charges for dataset: %s', dataset_name);
        return;
    end
end

CLogger.info('Building MGF charge cache for dataset: %s', dataset_name);
fid = obj.m_mapFid(dataset_name);
if fid == -1
    error('Failed to open MGF dataset: %s', dataset_name);
end

original_position = ftell(fid);
cleanup = onCleanup(@() fseek(fid, original_position, 'bof'));

status = fseek(fid, 0, 'bof');
if status ~= 0
    error('Failed to rewind MGF dataset: %s', dataset_name);
end

dataset_charge_map = containers.Map('KeyType', 'char', 'ValueType', 'double');
current_title = '';
current_charge = NaN;

while ~feof(fid)
    line = fgetl(fid);
    if ~ischar(line)
        continue;
    end

    if strncmp(line, 'BEGIN IONS', 10)
        current_title = '';
        current_charge = NaN;
    elseif startsWith(line, 'TITLE=')
        current_title = line(7:end);
    elseif startsWith(line, 'CHARGE=')
        current_charge = sscanf(line, 'CHARGE=%d+');
    elseif strcmp(line, 'END IONS')
        if ~isempty(current_title) && ~isempty(current_charge) && ~isnan(current_charge)
            dataset_charge_map(current_title) = current_charge;
            try
                scan_key = num2str(CMS2SpecNameUtils.parseMS2ScanNumber(current_title));
                if ~isequal(scan_key, current_title)
                    dataset_charge_map(scan_key) = current_charge;
                end
            catch
                % Keep raw TITLE mapping only when scan number is not parseable.
            end
        end
    end
end

cache_meta = struct( ...
    'source_name', mgf_info.name, ...
    'source_bytes', mgf_info.bytes, ...
    'source_datenum', mgf_info.datenum, ...
    'cache_version', 1);
save(cache_path, 'dataset_charge_map', 'cache_meta');
end


function cache_path = get_charge_cache_path(spec_dir_path, dataset_name)
% Store the charge cache next to PTMdecoder's MGF index cache.
[~, dataset_stem] = fileparts(char(dataset_name));
cache_path = fullfile(spec_dir_path, [dataset_stem, '_MGF_charge_map.mat']);
end


function is_valid = is_charge_cache_valid(cache_meta, mgf_info)
% Validate that the charge cache still matches its source MGF file.
is_valid = isstruct(cache_meta) ...
    && isfield(cache_meta, 'source_name') ...
    && isfield(cache_meta, 'source_bytes') ...
    && isfield(cache_meta, 'source_datenum') ...
    && isfield(cache_meta, 'cache_version') ...
    && strcmp(char(cache_meta.source_name), mgf_info.name) ...
    && isequal(cache_meta.source_bytes, mgf_info.bytes) ...
    && isequal(cache_meta.source_datenum, mgf_info.datenum) ...
    && isequal(cache_meta.cache_version, 1);
end
