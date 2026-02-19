function drawXIC(obj, dir_save, color_map, legend_map)
% Draw the XIC for gathered peptides using manually-checked rt range, to dir_save
% Input:
%   obj (CMSMSPepDeconv)
%       Processor instance
%   dir_save (1 x 1 char/string)
%       directory to save XIC figures
%   color_map (containers.Map or [])
%       color map (key: imp name, value: RGB 1x3)
%   legend_map (containers.Map or [])
%       legend map (key: imp name, value: display string)

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

report = CIMPQuantResultIO.read(checked_pep_path);
pep_rtrange_map = report.build_pep_rtrange_map();
if isempty(pep_rtrange_map)
    warning(['The checked peptide result file "', checked_pep_path, '" is empty!']);
end

if ~isfolder(dir_save)
    mkdir(dir_save);
end

%% Read the msms results and draw XICs
if isempty(obj.m_msms_res_path)
    each_PSM_results_path = fullfile(obj.m_outputDir, 'report_msms.txt');
else
    each_PSM_results_path = obj.m_msms_res_path;
end

% Indexing resources lazily
[obj, ~] = obj.ensureMs12DatasetIO();
[obj, mgf_created_here] = obj.ensureMgfDatasetIO();
if mgf_created_here
    cleanup_mgf = onCleanup(@() obj.m_cMgfDatasetIO.CloseAllFile());
end

% Initialize protein service lazily
[obj, ~] = obj.ensurePepProtService();

% Read and process
msms_reader = CMSMSResReader();
msms_result = msms_reader.read_from_msms_res_file(each_PSM_results_path);
print_progress = CPrintProgress(length(msms_result.Peptides));
pipeline = CIMPProcessingPipeline(obj.m_cMs12DatasetIO, obj.m_ms1_tolerance, obj.m_min_MSMS_num, obj.m_alpha, obj.m_resFilterThres);

fprintf('Drawing XIC...');
for idx_psf = 1:length(msms_result.Peptides)
    % Show progress
    print_progress = print_progress.update_show(idx_psf);
    rawIdentManager = obj.buildRawIdentManagerFromSpectrumList(msms_result.Peptides(idx_psf).spectrum_list);
    % Run gather
    pipeline.drawImpXicForBlock(rawIdentManager, pep_rtrange_map, dir_save, color_map, legend_map);
end

print_progress.last_update();
fprintf('done.\n');
end