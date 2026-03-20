function SetMap(obj)
% Create a dictionary mapping filenames to *_index, *_peaks
% build up spectral index for each dataset, including ms1_index, ms1_peaks, ms2_index
% Input:
%   obj (CMS12DatasetIO)
%       dataset IO instance

ms1_dataset_files=dir(fullfile(obj.m_strFoldname,'*.ms1'));
ms2_dataset_files=dir(fullfile(obj.m_strFoldname,'*.ms2'));

% begin building up
for i_ds=1:length(ms1_dataset_files)

    % get the MS1 info, and write to _MS1scan(peak)s.mat
    ms1_file_path = fullfile(ms1_dataset_files(i_ds).folder, ms1_dataset_files(i_ds).name);
    if 0==obj.load_MS1_file(ms1_file_path)   % Generate MS1 spectrum index file (_MS1scans.mat) and peak file (_MS1peaks.mat)
        error('Error in loading MS1 file: %s', ms1_file_path);
    end

    % get the MS2 info, and write to _MS2scan(peak)s.mat
    ms2_file_path = fullfile(ms2_dataset_files(i_ds).folder, ms2_dataset_files(i_ds).name);
    if 0==obj.load_MS1_MS2_mapping(ms2_file_path)   % Generate MS2 spectrum index file (_MS2scans.mat) and peak file (_MS2peaks.mat)
        error('Error in loading MS2 file: %s', ms2_file_path);
    end
    MS1_scanfile = fullfile(ms1_dataset_files(i_ds).folder,[ms1_dataset_files(i_ds).name(1:end-4),'_MS1scans.mat']);
    MS1_peakfile = fullfile(ms1_dataset_files(i_ds).folder,[ms1_dataset_files(i_ds).name(1:end-4),'_MS1peaks.mat']);
    MS2_scanfile = fullfile(ms2_dataset_files(i_ds).folder,[ms2_dataset_files(i_ds).name(1:end-4),'_MS2scans.mat']);
    % MS2_peakfile = fullfile(ms2_dataset_files(i_ds).folder,[ms2_dataset_files(i_ds).name(1:end-4),'_MS2peaks.mat']);
    load(MS1_scanfile, 'MS1_index');% MS1_index
    load(MS1_peakfile, 'MS1_peaks');% MS1_peaks
    load(MS2_scanfile, 'MS2_index');% MS2_index
    % load(MS2_peakfile, 'MS2_peaks');% MS2_peaks

    obj.m_mapNameMS1Index(ms1_dataset_files(i_ds).name(1:end-4))=MS1_index;
    obj.m_mapNameMS1Peaks(ms1_dataset_files(i_ds).name(1:end-4))=MS1_peaks;

    % Build sorted MS1 caches for fast m/z-window lookup.
    scan_boundaries = [0; MS1_index(:, 3); inf];
    peak_global_idx = (1:size(MS1_peaks, 1))';
    scan_idx_all = discretize(peak_global_idx, scan_boundaries);

    [sorted_mz, sort_order] = sort(MS1_peaks(:, 1), 'ascend');
    sorted_int = MS1_peaks(sort_order, 2);
    sorted_scan = scan_idx_all(sort_order);

    obj.m_mapNameMS1SortedMz(ms1_dataset_files(i_ds).name(1:end-4)) = sorted_mz;
    obj.m_mapNameMS1SortedInt(ms1_dataset_files(i_ds).name(1:end-4)) = sorted_int;
    obj.m_mapNameMS1SortedScan(ms1_dataset_files(i_ds).name(1:end-4)) = sorted_scan;

    obj.m_mapNameMS2Index(ms2_dataset_files(i_ds).name(1:end-4))=MS2_index;
    % obj.m_mapNameMS2Peaks(ms2_dataset_files(i_ds).name(1:end-4))=MS2_peaks;
end

end