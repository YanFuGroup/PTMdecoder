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

pep_res_reader = CPepResReader();
pep_res_reader = pep_res_reader.read_from_pep_res_file(checked_pep_path);
pep_rtrange_map = pep_res_reader.get_pep_rtrange_map();
if isempty(pep_rtrange_map)
    warning(['The checked peptide result file "', checked_pep_path, '" is empty!']);
end

mkdir(dir_save);

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
obj.CPepProtService = CPepProtService(obj.m_fastaFile, obj.m_regular_express, obj.m_filtered_res_file_path);

% Read and process
msms_reader = CMSMSResReader();
msms_result = msms_reader.read_from_msms_res_file(each_PSM_results_path);
print_progress = CPrintProgress(length(msms_result.Peptides));

fprintf('Drawing XIC...');
for idx_psf = 1:length(msms_result.Peptides)
    % Show progress
    print_progress = print_progress.update_show(idx_psf);
    % Get the peptide sequence
    peptide_sequence = msms_result.Peptides(idx_psf).peptide_sequence;
    % Get the protein name and position
    cell_prot_name_pos = obj.CPepProtService.get_protein_name_pos(peptide_sequence);
    % Initialize the CPepIsoGatherQuant object
    pepIsoGatherIMSLQ = CPepIsoGatherQuant(cell_prot_name_pos,obj.m_cMs12DatasetIO,...
        obj.m_resFilterThres,obj.m_ms1_tolerance,obj.m_alpha,'');
    % Get the spectrum list
    for idx_spec = 1:length(msms_result.Peptides(idx_psf).spectrum_list)
        % Get the dataset name and spectrum name
        dataset_name = msms_result.Peptides(idx_psf).spectrum_list(idx_spec).dataset_name;
        spectrum_name = msms_result.Peptides(idx_psf).spectrum_list(idx_spec).spectrum_name;
        peptidoform_strs = msms_result.Peptides(idx_psf).spectrum_list(idx_spec).peptidoform_list_str;
        peptidoform_abuns = msms_result.Peptides(idx_psf).spectrum_list(idx_spec).peptidoform_list_abun;
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