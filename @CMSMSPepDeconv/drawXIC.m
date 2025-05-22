function drawXIC(obj, dir_save, color_map, legend_map)
% Draw the XIC for gathered peptides using manually-checked rt range, to dir_save
% Input:
%   obj
%       the current object
%   dir_save
%       the directory to save the XIC figures
%   color_map
%       the color map
%   legend_map
%       the legend map

if nargin < 4
    legend_map = [];
end
if nargin < 3
    color_map = [];
end

%% Read the checked peptides and their XIC range
if isempty(obj.m_checked_peptides_res_path)
    checked_pep_path = fullfile(obj.m_outputDir, 'report_peptide_all_checked.txt');
else
    checked_pep_path = obj.m_checked_peptides_res_path;
end

fin = fopen(checked_pep_path, 'r');
if fin < 0
    error(['Cannot open the checked peptide level result:"',checked_pep_path,'"!']);
end
file_total_length = dir(checked_pep_path).bytes;
if file_total_length == 0
    fprintf(['Warning: The file "', checked_pep_path, '" is empty']);
end
print_progress = CPrintProgress(file_total_length);
fprintf('Reading the checked peptide results...');

% Initial the [mod_pep -> rt_range] data structure
pep_rtrange_map = containers.Map();

% Read the file and construct the [mod_pep -> rt_range] data structure
% Skip the first three lines (Header lines)
fgetl(fin);
fgetl(fin);
fgetl(fin);
key_mod_pep = '';
rt_ranges_temp = struct('rt_start',{},'rt_end',{},'check_label',{});
while ~feof(fin)
    strline = fgetl(fin);
    % Show progress
    now_bytes = ftell(fin);
    print_progress = print_progress.update_show(now_bytes);
    if isempty(strline)
        continue;
    elseif strline(1) == '@'
        % Record one retention time range line
        rt_ranges_temp = add_one_to_rt_ranges(rt_ranges_temp, strline);
    elseif strline(1) == '*'
        % Record one IMP line
        pep_rtrange_map = add_one_to_pep_rtrange_map(pep_rtrange_map,key_mod_pep,rt_ranges_temp);
        key_mod_pep = get_mod_pep_from_string(strline);
        rt_ranges_temp = struct('rt_start',{},'rt_end',{},'check_label',{});
    else
        % Record one protein-site line
        continue;
    end
end
% Record once more at the end of the file
pep_rtrange_map = add_one_to_pep_rtrange_map(pep_rtrange_map,key_mod_pep,rt_ranges_temp);
fclose(fin);
print_progress.last_update();
fprintf('done.\n');

%% Read the msms results and requantify the IMPs
if isempty(obj.m_msms_res_path)
    each_PSM_results_path = fullfile(obj.m_outputDir, 'report_msms.txt');
else
    each_PSM_results_path = obj.m_msms_res_path;
end

% fin = fopen(each_PSM_results_path, 'r');
% if fin < 0
%     error(['Cannot open the msms level result:"',each_PSM_results_path,'"!']);
% end

file_total_length = dir(each_PSM_results_path).bytes;
if file_total_length == 0
    fprintf(['Warning: The file "', each_PSM_results_path, '" is empty']);
end
% print_progress = CPrintProgress(file_total_length);

% Indexing the dataset IO
obj.m_cMs12DatasetIO = CMS12DatasetIO(obj.m_specPath,obj.m_ms1_tolerance);
obj.m_cMs12DatasetIO.SetMap();
obj.m_cMgfDatasetIO = CMgfDatasetIO;
obj.m_cMgfDatasetIO.Init(obj.m_specPath);
obj.m_cMgfDatasetIO.SetMap();
obj.m_cMgfDatasetIO.SetFidmap();

% Initial the fasta IO
obj.m_cFastaIO = CFastaIO(obj.m_fastaFile, obj.m_regular_express);

% Read and process
msms_reader = CMSMSResReader();
msms_reader = msms_reader.read_from_msms_res_file(each_PSM_results_path);
print_progress = CPrintProgress(length(msms_reader.m_peps_specs_forms));

fprintf('Re-quantifying at peptide level...');
for idx_psf = 1:length(msms_reader.m_peps_specs_forms)
    % Show progress
    print_progress = print_progress.update_show(idx_psf);
    % Get the peptide sequence
    peptide_sequence = msms_reader.m_peps_specs_forms(idx_psf).peptide_sequence;
    % Get the protein name and position
    cell_prot_name_pos = obj.m_cFastaIO.get_protein_name_pos(peptide_sequence);
    % Initialize the CPepIsoGatherQuant object
    pepIsoGatherIMSLQ = CPepIsoGatherQuant(cell_prot_name_pos,obj.m_cMs12DatasetIO,...
        obj.m_resFilterThres,obj.m_ms1_tolerance,obj.m_alpha,'');
    % Get the spectrum list
    for idx_spec = 1:length(msms_reader.m_peps_specs_forms(idx_psf).spectrum_list)
        % Get the dataset name and spectrum name
        dataset_name = msms_reader.m_peps_specs_forms(idx_psf).spectrum_list(idx_spec).dataset_name;
        spectrum_name = msms_reader.m_peps_specs_forms(idx_psf).spectrum_list(idx_spec).spectrum_name;
        peptidoform_strs = msms_reader.m_peps_specs_forms(idx_psf).spectrum_list(idx_spec).peptidoform_list_str;
        peptidoform_abuns = msms_reader.m_peps_specs_forms(idx_psf).spectrum_list(idx_spec).peptidoform_list_abun;
        % Get the profiles
        [isorts,c_ref_isointens,c_mz,cur_ch] = obj.getProfiles(dataset_name,spectrum_name);
        % Get the masses of IMPs
        lfMasses = get_masses_IMPs(peptidoform_strs,[obj.m_fixedModNameMass;obj.m_variableModNameMass]);
        % Append the quantification
        pepIsoGatherIMSLQ = pepIsoGatherIMSLQ.appendOneSpecQuant(dataset_name,isorts, ...
            c_ref_isointens,c_mz,cur_ch,peptidoform_strs,lfMasses,peptidoform_abuns);
    end
    % Run gather
    pepIsoGatherIMSLQ.drawGather(pep_rtrange_map, dir_save, color_map, legend_map);
end

print_progress.last_update();
fprintf('done.\n');
% fclose(fin);
obj.m_cMgfDatasetIO.CloseAllFile();
end



% Add one rt range to rt ranges
function rt_ranges = add_one_to_rt_ranges(rt_ranges, strline)
% Input:
%   rt_ranges
%       rt ranges of current IMP
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

% Add one rt range to the rt ranges struct
idx = length(rt_ranges) + 1;
rt_ranges(idx).rt_start = numbers(1);
rt_ranges(idx).rt_end = numbers(2);
rt_ranges(idx).check_label = numbers(4);
end



% Add one rt ranges to the pep_rtrange struct
function pep_rtrange_map = add_one_to_pep_rtrange_map(pep_rtrange_map,key_mod_pep,rt_ranges_temp)
% Input:
%   pep_rtrange_map
%       the modified peptide -> retention time range structure
%   key_mod_pep
%       the modified peptide name
%   rt_ranges_temp
%       the retention time ranges structure

% Skip the first line calling or empty-range peptide.
if isempty(key_mod_pep) || isempty(rt_ranges_temp)
    return;
end

% Create a new peptide record
pep_rtrange_map(key_mod_pep) = rt_ranges_temp;
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
function lfMasses = get_masses_IMPs(cstrPepIso,modNameMass)
% Input:
%   cstrPepIso
%       The IMP names
%   modNameMass
%       The cell of modification names, specificities and masses
% Output:
%   lfMasses
%       The masses of each IMPs
lfMasses = zeros(length(cstrPepIso),1);
for idx_iso = 1:length(cstrPepIso)
    % Split the sequence of peptide and the modification
    mod_seq = cstrPepIso{idx_iso};
    reg_exp = '\{(.*?)\}';
    [mod_str, seq_str] = regexp(mod_seq,reg_exp,'tokens','split');
    % Join the strings of sequence, delete the first and last "_" and count
    %   the masses.
    seq_str = strjoin(seq_str,'');
    seq_str([1,end]) = [];
    lfMasses(idx_iso) = sum(CConstant.vAAmass(seq_str-'A'+1));
    % Add the masses of modifications
    for idx_mod = 1:length(mod_str)
        is_notfound = true;
        for idx_mlist = 1:size(modNameMass,1)
            if isequal(modNameMass{idx_mlist,1},mod_str{idx_mod}{1})
                is_notfound = false;
                lfMasses(idx_iso) = lfMasses(idx_iso) + modNameMass{idx_mlist,3};
                break;
            end
        end
        if is_notfound
            error(['Unexpected modification is found: "',mod_str{idx_mod},'"!']);
        end
    end
    % Add the mass of water
    lfMasses(idx_iso) = lfMasses(idx_iso) + CConstant.hmass*2 + CConstant.omass;
end
end