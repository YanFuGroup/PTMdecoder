function obj = read_from_msms_res_file(obj, msms_res_path)
% Read from a msms result file
% Input:
%   msms_res_path
%       the path of the msms result file

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
        obj = obj.append_one_peptide(strline(3:end));
    elseif strline(1) == 'S'
        % Record one spectrum line
        strline = strsplit(strline, '\t');
        obj = obj.append_one_spectrum(strline{2}, strline{3});
    else
        % Record one peptidoform line
        strline = strsplit(strline, '\t');
        obj = obj.append_one_peptidoform(strline{1}, str2double(strline{2}));
    end
end

fclose(fin);

% Organize the data structure: delete the empty peptide, spectrum, and peptidoform
for i = length(obj.m_peps_specs_forms):-1:1
    if isempty(obj.m_peps_specs_forms(i).spectrum_list)
        obj.m_peps_specs_forms(i) = [];
    else
        for j = length(obj.m_peps_specs_forms(i).spectrum_list):-1:1
            if obj.m_peps_specs_forms(i).spectrum_list(j).peptidoform_num == 0
                obj.m_peps_specs_forms(i).spectrum_list(j) = [];
            else
                obj.m_peps_specs_forms(i).spectrum_list(j).peptidoform_list_str( ...
                    obj.m_peps_specs_forms(i).spectrum_list(j).peptidoform_num+1:end) = [];
                obj.m_peps_specs_forms(i).spectrum_list(j).peptidoform_list_abun( ...
                    obj.m_peps_specs_forms(i).spectrum_list(j).peptidoform_num+1:end) = [];
            end
        end
        if isempty(obj.m_peps_specs_forms(i).spectrum_list)
            obj.m_peps_specs_forms(i) = [];
        end
    end
end

end

