function resultObj = read_from_msms_res_file(~, msms_res_path)
% Read from a msms result file
% Input:
%   msms_res_path
%       the path of the msms result file
% Output:
%   resultObj
%       CMSMSResult object containing the parsed data

% Initialize result object
resultObj = CMSMSResult();

fin = fopen(msms_res_path, 'r');
if fin < 0
    error(['Cannot open the msms level result:"',msms_res_path,'"!']);
end
file_total_length = dir(msms_res_path).bytes;
if file_total_length == 0
    error(['Warning: The file "',msms_res_path,'" is empty!']);
    return; %#ok<UNRCH> 
end

% Read the file
while(~feof(fin))
    strline = fgetl(fin);
    if isempty(strline)
        continue;
    end
    if strline(1) == 'P'
        % Record one peptide line
        resultObj.addPeptide(strline(3:end));
    elseif strline(1) == 'S'
        % Record one spectrum line
        strline = strsplit(strline, '\t');
        resultObj.addSpectrum(strline{2}, strline{3});
    else
        % Record one peptidoform line
        strline = strsplit(strline, '\t');
        resultObj.addPeptidoform(strline{1}, str2double(strline{2}));
    end
end

fclose(fin);

% Organize the data structure: compress buffers
resultObj.compress();

end

