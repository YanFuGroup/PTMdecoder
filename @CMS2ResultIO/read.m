function resultObj = read(msms_res_path)
% Read from a msms result file
% Input:
%   msms_res_path (1 x 1 char/string)
%       path of the msms result file
% Output:
%   resultObj (CMS2Result)
%       CMS2Result object containing the parsed data

% Initialize result object
resultObj = CMS2Result();

fin = fopen(msms_res_path, 'r');
if fin < 0
    error(['Cannot open the msms level result:"', msms_res_path, '"!']);
end
cleanup = onCleanup(@() fclose(fin));

file_total_length = dir(msms_res_path).bytes;
if file_total_length == 0
    error(['Warning: The file "', msms_res_path, '" is empty!']);
end

% Read the file
while(~feof(fin))
    strline = fgetl(fin);
    if isempty(strline)
        continue;
    end
    tokens = strsplit(strline, '\t');
    tag = tokens{1};
    if strcmp(tag, 'P')
        % Record one peptide line (allow repeated P lines for the same sequence)
        if numel(tokens) < 2
            error('CMS2ResultIO:InvalidPeptideLine', 'Invalid peptide line: %s', strline);
        end
        resultObj.addOrSelectPeptide(tokens{2});
    elseif strcmp(tag, 'S')
        % Record one spectrum line
        if numel(tokens) < 3
            error('CMS2ResultIO:InvalidSpectrumLine', 'Invalid spectrum line: %s', strline);
        end
        resultObj.addSpectrum(tokens{2}, tokens{3});
    else
        % Record one peptidoform line
        if numel(tokens) < 2
            error('CMS2ResultIO:InvalidPeptidoformLine', 'Invalid peptidoform line: %s', strline);
        end
        resultObj.addPeptidoform(tokens{1}, str2double(tokens{2}));
    end
end

% Organize the data structure: compress buffers
resultObj.compress();
end
