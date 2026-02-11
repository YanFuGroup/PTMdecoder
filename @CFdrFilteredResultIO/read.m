function result = read(filePath)
% Read FDR filtered result table from file.
% Input:
%   filePath (1 x 1 char/string)
%       Path to the FDR filtered result table.
% Output:
%   result (CFdrFilteredResult)
%       Container holding parsed entries.
%
% The reader expects a tab-delimited table with a header line followed by
% data rows in the Mascot-style column order:
%   Site, DatasetName, Scan, Spectrum, Charge, Calc_neutral_pepmass,
%   precursor_neutral_mass, massdiff, num_match_ions, peptide, protein,
%   modification, modificationlocation, Score.

fin = fopen(filePath, 'r');
if fin == -1
    error('Cannot open the FDR filtered result file: "%s"!', filePath);
end
cleanup = onCleanup(@() fclose(fin));

fileInfo = dir(filePath);
progress_printer = CPrintProgress(fileInfo.bytes);

header = fgetl(fin);
if ~ischar(header) && ~isstring(header)
    error('Empty FDR filtered result file: "%s"!', filePath);
end

entries = CFdrFilteredResult.emptyEntries();
lineNo = 1;
while ~feof(fin)
    strline = fgetl(fin);
    lineNo = lineNo + 1;
    progress_printer = progress_printer.update_show(ftell(fin));

    if isempty(strline)
        continue;
    end

    segment = regexp(strline, '\t| +', 'split');
    segment = segment(~cellfun('isempty', segment));

    if numel(segment) < 14
        error('Invalid FDR filtered line %d in "%s": %s', lineNo, filePath, strline);
    end

    entry = struct( ...
        'Site', segment{1}, ...
        'DatasetName', segment{2}, ...
        'Scan', segment{3}, ...
        'Spectrum', segment{4}, ...
        'Charge', segment{5}, ...
        'Calc_neutral_pepmass', segment{6}, ...
        'precursor_neutral_mass', segment{7}, ...
        'massdiff', segment{8}, ...
        'num_match_ions', segment{9}, ...
        'peptide', segment{10}, ...
        'protein', segment{11}, ...
        'modification', segment{12}, ...
        'modificationlocation', segment{13}, ...
        'Score', segment{14} ...
    );
    entries(end + 1) = entry; %#ok<AGROW>
end

progress_printer.last_update();
result = CFdrFilteredResult(entries);
end
