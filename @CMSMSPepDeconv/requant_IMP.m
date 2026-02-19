function obj = requant_IMP(obj)
% Re-quantify the IMPs using checked XIC peaks
% Input:
%   obj (CMSMSPepDeconv)
%       Processor instance
% Output:
%   obj (CMSMSPepDeconv)
%       Updated instance

%% Read the checked peptides and their XIC range
if isempty(obj.m_checked_peptides_res_path)
    checked_pep_path = fullfile(obj.m_outputDir, 'report_peptide_all_checked.txt');
else
    checked_pep_path = obj.m_checked_peptides_res_path;
end

report = CIMPQuantResultIO.read(checked_pep_path);
pep_rtrange_map = report.build_pep_rtrange_map();
if isempty(pep_rtrange_map)
    error(['The checked peptide result file "', checked_pep_path, '" is empty!']);
end

if ~isfolder(obj.m_outputDir)
    mkdir(obj.m_outputDir);
end

%% Read the msms results and requantify the IMPs
if isempty(obj.m_msms_res_path)
    each_PSM_results_path = fullfile(obj.m_outputDir, 'report_msms.txt');
else
    each_PSM_results_path = obj.m_msms_res_path;
end

output_path = fullfile(obj.m_outputDir, 'report_peptide_all_requant.txt');
report = CIMPQuantReport();

% Indexing resources lazily
[obj, ~] = obj.ensureMs12DatasetIO();
[obj, mgf_created_here] = obj.ensureMgfDatasetIO();
if mgf_created_here
    cleanup_mgf = onCleanup(@() obj.m_cMgfDatasetIO.CloseAllFile());
end

% check and create a new output file

% Initialize protein service lazily
[obj, ~] = obj.ensurePepProtService();

% Read and process
msms_reader = CMSMSResReader();
msms_result = msms_reader.read_from_msms_res_file(each_PSM_results_path);
print_progress = CPrintProgress(length(msms_result.Peptides));
if ~isempty(obj.m_taskParam)
    overrides = struct('ms1_tolerance', obj.m_ms1_tolerance, 'minMSMSnum', obj.m_min_MSMS_num, ...
        'alpha', obj.m_alpha, 'resFilterThres', obj.m_resFilterThres);
    pipeline_cfg = CIMPProcessingPipelineConfig.fromTaskParam(obj.m_taskParam, obj.m_cMs12DatasetIO, overrides);
else
    pipeline_cfg = CIMPProcessingPipelineConfig(obj.m_cMs12DatasetIO, obj.m_ms1_tolerance, obj.m_min_MSMS_num, obj.m_alpha, obj.m_resFilterThres);
end
pipeline = CIMPProcessingPipeline(pipeline_cfg);

fprintf('Re-quantifying at peptide level...')
for idx_psf = 1:length(msms_result.Peptides)
    % Show progress
    print_progress = print_progress.update_show(idx_psf);
    % Get the peptide sequence
    peptide_sequence = msms_result.Peptides(idx_psf).peptide_sequence;
    % Get the protein name and position
    cell_prot_name_pos = obj.CPepProtService.get_protein_name_pos(peptide_sequence);
    rawIdentManager = obj.buildRawIdentManagerFromSpectrumList(msms_result.Peptides(idx_psf).spectrum_list);
    % Run gather
    block = pipeline.requantifyPeptideBlock(cell_prot_name_pos, rawIdentManager, pep_rtrange_map);
    report = report.append_block(block);

end

print_progress.last_update();
fprintf('done.\n');
CIMPQuantResultIO.write(report, output_path);
end