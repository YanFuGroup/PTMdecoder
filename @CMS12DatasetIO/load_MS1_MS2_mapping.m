function success = load_MS1_MS2_mapping(~, ms2_fullfile)
% Parse MS2 file to get mapping between MS2 scan and its precursor MS1 scan
% Output file: *_MS2scans.mat containing MS2_index variable
% MS2_index structure:
%   Column 1: MS1 scan number (PrecursorScan)
%   Column 2: MS2 scan number (from S line)
% Input:
%   ms2_fullfile (1 x 1 char/string)
%       MS2 file path
% Output:
%   success (1 x 1 double/int)
%       1 for success, 0 for failure

success = 0;

% Check input file
if ~exist(ms2_fullfile, 'file')
    fprintf('%s: does not exist!\n', ms2_fullfile);
    return;
end

% Prepare output file path
[datapath, dataname] = fileparts(ms2_fullfile);
MS2_scanfile = fullfile(datapath, [dataname, '_MS2scans.mat']);

% Check if output file already exists
if exist(MS2_scanfile, 'file')
    success = 1;
    return;
end

% Open file
fid = fopen(ms2_fullfile, 'r');
if fid == -1
    fprintf('Cannot open: %s\n', ms2_fullfile);
    return;
end

fprintf('Processing MS2 mapping: %s\n', dataname);

% Initialize variables
% Estimate number of scans to pre-allocate. 1.5e5 is used in original code
maxMS2num = 1.5e5;
MS2_index = zeros(maxMS2num, 2); 
scan_count = 0;

% Keywords
key_S = 'S';
key_Precursor = 'I	PrecursorScan';

len_S = length(key_S);
len_Precursor = length(key_Precursor);

curr_ms2_scan = 0;
curr_ms1_scan = 0;
found_S = false;
found_P = false;

tline = fgets(fid);
while ischar(tline)
    % Check for 'S' line: S  scan  scan  mz
    if strncmp(tline, key_S, len_S)
        % If we found a previous pair, save it
        if found_S && found_P
            scan_count = scan_count + 1;
            MS2_index(scan_count, :) = [curr_ms1_scan, curr_ms2_scan];
        end
        
        % Reset flags for new spectrum
        found_S = true;
        found_P = false;
        
        % Parse S line
        % Format: S [tab] ScanNum [tab] ScanNum [tab] PrecursorMz
        parts = textscan(tline, '%s %f %f %f');
        if ~isempty(parts{2})
            curr_ms2_scan = parts{2};
        end
        
    elseif strncmp(tline, key_Precursor, len_Precursor)
        % Check for 'I PrecursorScan' line
        % Format: I [tab] PrecursorScan [tab] ScanNum
        % Using sscanf to extract the number at the end
        % The line is "I	PrecursorScan	2"
        
        % Extract the scan number.
        % We know the prefix length.
        str_val = tline(len_Precursor+2:end); 
        curr_ms1_scan = str2double(str_val);
        
        if ~isnan(curr_ms1_scan)
            found_P = true;
        end
    end
    
    tline = fgets(fid);
end

% Handle the last entry
if found_S && found_P
    scan_count = scan_count + 1;
    MS2_index(scan_count, :) = [curr_ms1_scan, curr_ms2_scan];
end

fclose(fid);

% Trim the result
MS2_index = MS2_index(1:scan_count, :);

% Save to .mat file
save(MS2_scanfile, 'MS2_index');

fprintf('Finished. Found %d scans.\n', scan_count);

success = 1;
end
