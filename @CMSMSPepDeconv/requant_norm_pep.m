function obj = requant_norm_pep(obj)
% Re-quantify the normalization peptides using checked XIC peaks

%% Read the checked peptides and their XIC range
if isempty(obj.m_checked_peptides_res_path)
    checked_pep_path = fullfile(obj.m_outputDir, 'peptide4normalization_checked.txt');
else
    checked_pep_path = obj.m_checked_peptides_res_path;
end

fin = fopen(checked_pep_path, 'r');
if fin < 0
    error(['Cannot open the checked normalization peptides result:"',checked_pep_path,'"!']);
end

file_total_length = dir(checked_pep_path).bytes;
if file_total_length == 0
    fprintf(['Warning: The file "', checked_pep_path, '" is empty\n']);
end
print_progress = CPrintProgress(file_total_length);
fprintf('Reading the checked peptide results...');

% Initial the [mod_pep -> rt_range] data structure
pep_rtrange_map = containers.Map();
pep_prot_map = containers.Map();

% Read the file and construct the [mod_pep -> rt_range] data structure
% Skip the first three lines (Header lines)
% TODO: Use CPepResReader
fgetl(fin);
fgetl(fin);
fgetl(fin);
prot_name = '';
while ~feof(fin)
    strline = fgetl(fin);
    % Show progress
    print_progress = print_progress.update_show(ftell(fin));
    if isempty(strline)
        continue;
    elseif strline(1) == '@'
        % Record one retention time range line
        [rt_left, rt_right, check_label] = get_rt_range_check_label(strline);
        pep_rtrange_map(key_mod_pep) = [pep_rtrange_map(key_mod_pep); struct('rt_start',rt_left,'rt_end',rt_right,'check_label',check_label)];
    elseif strline(1) == '*'
        % Record one IMP line
        key_mod_pep = get_mod_pep_from_string(strline);
        pep_prot_map(key_mod_pep) = prot_name;
        if ~isempty(key_mod_pep)
            pep_rtrange_map(key_mod_pep) = struct('rt_start',{},'rt_end',{},'check_label',{});
        end
    else
        % Record one protein-site line
        % Only retain the name of first protein
        segment = regexp(strline,',','split');
        prot_name = segment{1};
        continue;
    end
end
% Record once more at the end of the file
fclose(fin);
print_progress.last_update();
fprintf('done.\n');

mkdir(obj.m_outputDir);

%% Requantify the normalization peptides
output_path = fullfile(obj.m_outputDir, 'peptide4normalization_requant.txt');
fout = fopen(output_path, 'w');
if fout < 0
    error(['Cannot open the re-quantification result file:"',output_path,'"!']);
end
fprintf(fout,'Protein_name,Peptide_start_position_on_protein;\n');
fprintf(fout,'*\tIMP\tCharge\tDataset\tMass_center\tLow_mass_bound\tHigh_mass_bound\tPeak_area\n');
fprintf(fout,'@\tRT_start\tRT_end\tProportion\tCheck_label\n');
fclose(fout);

% Indexing the dataset IO
obj.m_cMs12DatasetIO = CMS12DatasetIO(obj.m_specPath,obj.m_ms1_tolerance);
obj.m_cMs12DatasetIO.SetMap();
obj.m_cMgfDatasetIO = CMgfDatasetIO;
obj.m_cMgfDatasetIO.Init(obj.m_specPath);
obj.m_cMgfDatasetIO.SetMap();
obj.m_cMgfDatasetIO.SetFidmap();

% Read and process
fin = fopen(checked_pep_path, 'r');
if fin < 0
    error(['Cannot open the checked normalization peptides result:"',checked_pep_path,'"!']);
end

file_total_length = dir(checked_pep_path).bytes;
if file_total_length == 0
    fprintf(['Warning: The file "', checked_pep_path, '" is empty\n']);
end
print_progress = CPrintProgress(file_total_length);
fprintf('Re-quantifying at peptide level...')

% Read the file and construct the [mod_pep -> rt_range] data structure
% Skip the first three lines (Header lines)
% TODO: Use CPepResReader
fgetl(fin);
fgetl(fin);
fgetl(fin);
is_ready = false;
while ~feof(fin)
    strline = fgetl(fin);
    % Show progress
    print_progress = print_progress.update_show(ftell(fin));
    if isempty(strline)
        continue;
    elseif strline(1) == '@'
        % Record one retention time line, complete a psm result
        rt_median = get_median_rt(strline);
        pepIsoGatherIMSLQ = pepIsoGatherIMSLQ.appendOneSpecQuant(...
            mgf_name, rt_median, 1,lfMz, current_charge, ...
            {current_peptide}, lfMass, 1);
    elseif strline(1) == '*'
        % Record one peptide line
        if is_ready
            pepIsoGatherIMSLQ.rerunGather_quant(pep_rtrange_map);
        end
        [mgf_name, current_charge, current_peptide] = ...
            get_information_from_peptide_line(strline);
        lfMass = get_mass_peptide(current_peptide);
        lfMz = (lfMass+current_charge*CConstant.pmass)/current_charge;
        current_key = [current_peptide,'_+',num2str(current_charge),'_',mgf_name];
        pepIsoGatherIMSLQ = CPepIsoGatherQuant({pep_prot_map(current_key),-1}, ...
            obj.m_cMs12DatasetIO,obj.m_resFilterThres,obj.m_ms1_tolerance, ...
            obj.m_alpha,output_path);
        is_ready = true;
    else
        % Record one protein-site line
        if is_ready
            pepIsoGatherIMSLQ.rerunGather_quant(pep_rtrange_map);
            is_ready = false;
        end
    end
end
% Record once more at the end of the file
pepIsoGatherIMSLQ.rerunGather_quant(pep_rtrange_map);
fclose(fin);
print_progress.last_update();
fprintf('done.\n');

obj.m_cMgfDatasetIO.CloseAllFile();
end



%% Other functions

% Get the median of rt bound
function rt = get_median_rt(strline)
% Input:
%   strline
%       the current line

reg_exp_pat = '\d*\.?\d+';
% Use regexp to find all matches
matches = regexp(strline, reg_exp_pat, 'match');
numbers = str2double(matches);

% Check the size of the input numbers
if length(numbers) < 4
    error(['The line: "',strline,'" representing the RT ranges are in an unexpected format!']);
end

% Add one rt to identification result
rt = (numbers(1) + numbers(2)) / 2;
end



% Get some information for this peptide
function [mgf_name, current_charge, current_peptide] = ...
    get_information_from_peptide_line(strline)
% Input:
%   strline
%       One peptide line
% Output:
%   mgf_name
%       Name of mgf
%   current_charge
%       Charge of this PSM
%   current_peptide
%       Peptide sequence of this PSM
%   peptide_line_minus_1
%       String for output in this peptide line minus the last element

segment = regexp(strline,'\t','split');
% Need the modified peptide (2), charge (3) and dataset name (4)
mgf_name = segment{4};
current_charge = str2double(segment{3}(2:end)); % only use the number
current_peptide = segment{2};
end



% Get the rt range and check label
function [rt_left, rt_right, check_label] = get_rt_range_check_label(strline)
% Input:
%   strline
%       the current line
% Output:
%   rt_left
%       the left bound of rt range
%   rt_right
%       the right bound of rt range
%   check_label
%       the check label

reg_exp_pat = '\d*\.?\d+';
% Use regexp to find all matches
matches = regexp(strline, reg_exp_pat, 'match');
numbers = str2double(matches);

% Check the size of the input numbers
if length(numbers) < 4
    error(['The line: "',strline,'" representing the RT ranges are in an unexpected format!']);
end

% Get the rt range and check label
rt_left = numbers(1);
rt_right = numbers(2);
check_label = numbers(4);
end



% Get the key of mod peptide with a string
function key_mod_pep = get_mod_pep_from_string(strline)
% Input:
%   strline
%       the input string, read from report_peptide_all_checked
% Output:
%   key_mod_pep
%       the modified peptide strings

segment = regexp(strline,'\t','split');
% Need the modified peptide (2), charge (3) and dataset name (4)
key_mod_pep = [segment{2},'_',segment{3},'_',segment{4}];
end



% Get the mass of each IMPs
function lfMass = get_mass_peptide(pep_seq)
% Input:
%   pep_seq
%       The peptide sequence
% Output:
%   lfMass
%       The mass of the peptide

% Add the mass of each amino acid
lfMass = sum(CConstant.vAAmass(pep_seq-'A'+1));

% Add the mass of water
lfMass = lfMass + CConstant.hmass*2 + CConstant.omass;
end