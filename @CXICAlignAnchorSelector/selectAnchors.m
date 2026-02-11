function anchors = selectAnchors(~, fdr_file_path, ms12DatasetIO, options)
% Select unmodified peptide anchors from FDR filtered results.
% Input:
%   fdr_file_path (1 x 1 char/string)
%       FDR filtered result file path
%   ms12DatasetIO (CMS12DatasetIO)
%       MS1/MS2 dataset IO for RT lookup
%   options (struct, optional)
%       Selector options (min_psm)
% Output:
%   anchors (struct array)
%       Fields: peptide, raw_name, rt, count

if nargin < 4
    options = struct();
end
min_psm = COptionUtils.get(options, 'min_psm', 3);

fdr_results = CFdrFilteredResultIO.read(fdr_file_path);
entries = fdr_results.entries;

anchor_map = containers.Map('KeyType', 'char', 'ValueType', 'any');

for idx = 1:length(entries)
    entry = entries(idx);
    if ~isequal(entry.modification, '-')
        continue;
    end
    rt = resolve_entry_rt(entry, ms12DatasetIO);
    if isnan(rt)
        continue;
    end
    key = [entry.peptide, '|', entry.DatasetName];
    if ~isKey(anchor_map, key)
        anchor_map(key) = rt;
    else
        anchor_map(key) = [anchor_map(key), rt];
    end
end

keys = anchor_map.keys;
anchors = struct('peptide', {}, 'raw_name', {}, 'rt', {}, 'count', {});
for idx = 1:numel(keys)
    key = keys{idx};
    parts = regexp(key, '\|', 'split');
    peptide = parts{1};
    raw_name = parts{2};
    rts = anchor_map(key);
    if numel(rts) < min_psm
        continue;
    end
    anchors(end+1) = struct( ...
        'peptide', peptide, ...
        'raw_name', raw_name, ...
        'rt', median(rts), ...
        'count', numel(rts)); %#ok<AGROW>
end
end

function rt = resolve_entry_rt(entry, ms12DatasetIO)
% Resolve RT for one FDR entry from MS2->MS1 mapping.
% Input:
%   entry (struct)
%       FDR entry with DatasetName and Scan
%   ms12DatasetIO (CMS12DatasetIO)
%       MS1/MS2 dataset IO
% Output:
%   rt (double)
%       Retention time in minutes, NaN if missing

rt = NaN;
try
    mgf_name = entry.DatasetName;
    mgf_stem = erase(mgf_name, '.mgf');
    ms2_name = ms12DatasetIO.m_cMsFileMapper.get_ms2_stem(mgf_stem);
    MS2_index = ms12DatasetIO.m_mapNameMS2Index(ms2_name);
    curr_MS2_scan = str2double(entry.Scan);
    tmp_idx = MS2_index(:,2) == curr_MS2_scan;
    if ~any(tmp_idx)
        return;
    end
    MS1Scan = MS2_index(tmp_idx, 1);
    MS1_index = ms12DatasetIO.m_mapNameMS1Index(ms2_name);
    ino = find(MS1_index(:,1) == MS1Scan, 1);
    if isempty(ino)
        return;
    end
    rt = MS1_index(ino, 2);
catch
    rt = NaN;
end
end

