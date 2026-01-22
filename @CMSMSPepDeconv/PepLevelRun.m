function obj = PepLevelRun(obj)
% Run the quantification at peptide level
% Read the results of PSM level and quantify at peptide level
% Input:
%   obj - The CMSMSPepDeconv object
% Output:
%   obj - The CMSMSPepDeconv object

obj.check_whether_ms12_mgf_name_match();

% Load protein sequences from fasta file
if isempty(obj.m_cFastaIO)
    obj.m_cFastaIO = CFastaIO(obj.m_fastaFile, obj.m_regular_express, obj.m_filtered_res_file_path);
end

% Indexing the mgf
if isempty(obj.m_cMgfDatasetIO)
    obj.m_cMgfDatasetIO = CMgfDatasetIO;
    obj.m_cMgfDatasetIO.Init(obj.m_specPath);
    obj.m_cMgfDatasetIO.SetMap();
    obj.m_cMgfDatasetIO.SetFidmap();
    need_release_mgf_index = true;
else
    need_release_mgf_index = false;
end

% Indexing the ms1/ms2
if isempty(obj.m_cMs12DatasetIO)
    obj.m_cMs12DatasetIO = CMS12DatasetIO(obj.m_specPath,obj.m_ms1_tolerance);
    obj.m_cMs12DatasetIO.SetMap();
end

% Check the report_msms.txt file
if isempty(obj.m_msms_res_path)
    each_PSM_results_path = fullfile(obj.m_outputDir, 'report_msms.txt');
else
    each_PSM_results_path = obj.m_msms_res_path;
end
% fin=fopen(each_PSM_results_path,'r');
% if 0 >= fin
%     error(['Can not open file: ',each_PSM_results_path]);
% end
file_total_length = dir(each_PSM_results_path).bytes;
if file_total_length <= 0
    error(['The file "', each_PSM_results_path,'" is empty.'])
end
% print_progress = CPrintProgress(file_total_length);

% Check and create a new output file
each_peptide_results_path = fullfile(obj.m_outputDir,'report_peptide_all.txt');
fout = fopen(each_peptide_results_path,'w');
if fout <= 0
    error(['Cannot open the peptide level report file ',each_peptide_results_path]);
end
fprintf(fout,'Protein_name,Peptide_start_position_on_protein;\n');
fprintf(fout,'*\tIMP\tCharge\tDataset\tMass_center\tLow_mass_bound\tHigh_mass_bound\tPeak_area\n');
fprintf(fout,'@\tRT_start\tRT_end\tProportion\tCheck_label\n');
fclose(fout);

% Read and process
msms_reader = CMSMSResReader();
msms_reader = msms_reader.read_from_msms_res_file(each_PSM_results_path);
print_progress = CPrintProgress(length(msms_reader.m_peps_specs_forms));

fprintf('Quantifying at peptide level...')
for idx_psf = 1:length(msms_reader.m_peps_specs_forms)
    % Show progress
    print_progress = print_progress.update_show(idx_psf);
    % Get the peptide sequence
    peptide_sequence = msms_reader.m_peps_specs_forms(idx_psf).peptide_sequence;
    % Get the protein name and position
    cell_prot_name_pos = obj.m_cFastaIO.get_protein_name_pos(peptide_sequence);
    % Initialize the CPepIsoGatherQuant object
    pepIsoGatherIMSLQ = CPepIsoGatherQuant(cell_prot_name_pos,obj.m_cMs12DatasetIO,...
        obj.m_resFilterThres,obj.m_ms1_tolerance,obj.m_alpha,each_peptide_results_path,obj.m_min_MSMS_num);
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
    pepIsoGatherIMSLQ.runGather();
end
print_progress.last_update();
fprintf('done.\n');

if need_release_mgf_index
    obj.m_cMgfDatasetIO.CloseAllFile();
end
end



%% Other functions

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