function obj = read_from_pep_res_file(obj, pep_res_path)
% Read peptide results from a report file
% Input:
%   obj (CPepResReader)
%       Peptide result reader instance
%   pep_res_path (1 x 1 char/string)
%       Path to the peptide result file
% Output:
%   obj (CPepResReader)
%       Updated reader instance
fin = fopen(pep_res_path, 'r');
if fin < 0
    error(['Cannot open the checked peptide level result:"',pep_res_path,'"!']);
end

obj.m_pep_res = containers.Map();  % Initialize the overall peptide results
obj.m_prot_pep_res = {};  % Initialize the protein-peptide mapping as empty cell array

% Read the file and construct the [mod_pep -> rt_range] data structure
% Skip the first three lines (Header lines)
fgetl(fin);
fgetl(fin);
fgetl(fin);

% Initialize variables
key_mod_pep = '';
peptidoform_name = '';
charge_state = 0;
dataset_name = '';
mean_mz = 0;
lb_mz = 0;
ub_mz = 0;
quant_value = 0;
rt_ranges_temp = struct('rt_start',{},'rt_end',{},'ratio',{},'check_label',{});
ori_pep_line = '';

while ~feof(fin)
    strline = fgetl(fin);
    if isempty(strline)
        continue;
    elseif strline(1) == '@'
        % Record one retention time range line
        [rt_start, rt_end, ratio, check_label] = extract_rt_ranges(strline);
        rt_ranges_temp(end+1) = struct('rt_start', rt_start, 'rt_end', rt_end, 'ratio', ratio, 'check_label', check_label); %#ok<AGROW> 
    elseif strline(1) == '*'
        % Record one IMP line
        obj = obj.append_one_pep(key_mod_pep, peptidoform_name, charge_state, dataset_name, ...
            mean_mz, lb_mz, ub_mz, quant_value, rt_ranges_temp, ori_pep_line);
        [key_mod_pep, peptidoform_name, charge_state, dataset_name, mean_mz, lb_mz, ub_mz, quant_value] = get_pep_info_from_line(strline);
        % Add peptide to the most recent protein's peptide list
        if ~isempty(obj.m_prot_pep_res)
            obj.m_prot_pep_res{end, 2}{end+1} = peptidoform_name;
        end
        ori_pep_line = strline;  % Store the original peptide line
        rt_ranges_temp = struct('rt_start',{},'rt_end',{},'ratio',{},'check_label',{});
    else
        % Record one protein-site line
        obj.m_prot_pep_res{end+1, 1} = strline;  % Store protein name
        obj.m_prot_pep_res{end, 2} = {};           % Initialize peptide list
    end
end
% Record once more at the end of the file
obj = obj.append_one_pep(key_mod_pep, peptidoform_name, charge_state, dataset_name, ...
    mean_mz, lb_mz, ub_mz, quant_value, rt_ranges_temp, ori_pep_line);
fclose(fin);

end



% Add one rt range to rt ranges
function [rt_start, rt_end, ratio, check_label] = extract_rt_ranges(strline)
% Input:
%   strline (1 x 1 char/string)
%       the current line
% Output:
%   rt_start (1 x 1 double) minutes
%   rt_end (1 x 1 double) minutes
%   ratio (1 x 1 double)
%   check_label (1 x 1 double/int)

reg_exp_pat = '\d*\.?\d+';
% Use regexp to find all matches
matches = regexp(strline, reg_exp_pat, 'match');
numbers = str2double(matches);

% Check the size of the input numbers
if length(numbers) < 4
    error(['The line: "',strline,'" representing the RT ranges are in an unexpected format!']);
end

% Extract RT ranges from the string
rt_start = numbers(1);
rt_end = numbers(2);
ratio = numbers(3);
check_label = numbers(4);
end



% Get the key of mod peptide with a string
function [key_mod_pep, peptidoform_name, charge_state, dataset_name, ...
    mean_mz, lb_mz, ub_mz, quant_value] = get_pep_info_from_line(strline)
% Input:
%   strline (1 x 1 char/string)
%       the input string, read from report_peptide_all_checked
% Output:
%   key_mod_pep (1 x 1 char/string)
%       modified peptide key (peptide_charge_dataset)
%   peptidoform_name (1 x 1 char/string)
%       modified peptide name
%   charge_state (1 x 1 double/int)
%       charge state of the peptide
%   dataset_name (1 x 1 char/string)
%       dataset name
%   mean_mz (1 x 1 double) m/z
%       mean m/z value
%   lb_mz (1 x 1 double) m/z
%       lower bound m/z value
%   ub_mz (1 x 1 double) m/z
%       upper bound m/z value
%   quant_value (1 x 1 double)
%       quantification value of the peptide

segment = regexp(strline,'\t','split');
% Key is the combination of the modified peptide (2), charge_state (3) and dataset name (4)
key_mod_pep = [segment{2},'_',segment{3},'_',segment{4}];
peptidoform_name = segment{2};
charge_state = segment{3};
dataset_name = segment{4};
mean_mz = str2double(segment{5});
lb_mz = str2double(segment{6});
ub_mz = str2double(segment{7});
quant_value = str2double(segment{8});
end