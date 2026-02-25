function obj = runImpAlignRequantLevel(obj, align_strategy, align_options)
% Quantify -> align XIC -> requantify IMPs
% Input:
%   obj (CMSMSPepDeconv)
%       Processor instance
%   align_strategy (CRunAlignStrategy)
%       Alignment strategy (reference-run or pairwise)
%   align_options (struct, optional)
%       Alignment options (paths and parameters):
%         - msms_res_path (char/string)
%             Path to MSMS result report. Default: obj.m_msms_res_path or fullfile(obj.m_outputDir, 'report_msms.txt').
%         - alignment_report_path (char/string)
%             Output path for alignment report. Default: fullfile(obj.m_outputDir, 'report_alignment.txt').
%         - requant_output_path (char/string)
%             Output path for re-quant report. Default: fullfile(obj.m_outputDir, 'report_peptide_all_requant_aligned.txt').
%         - min_psm (1 x 1 double)
%             Minimum PSM count per peptide to form anchors. Default: 3.
%         - num_bins (1 x 1 double)
%             Number of bins for local alignment offsets. Default: 5.
%         - min_per_bin (1 x 1 double)
%             Minimum anchors per bin for local offset. Default: 5.
%         - outlier_k (1 x 1 double)
%             Outlier removal k (MAD/STD) for transform fitting. Default: 3.
%         - outlier_method (char/string)
%             Outlier method ('mad' or 'std'). Default: 'mad'.
%         - rt_sigma (1 x 1 double)
%             RT Gaussian sigma (minutes) for peak selection. Default: 0.5.
%         - max_rt_residual (1 x 1 double)
%             Max allowed RT residual (minutes) for peak pairing. Default: model.rt_residual_threshold.
%         - dead_time_min (1 x 1 double)
%             Min allowed RT start (minutes) for peaks. Default: [] (disabled).
% Output:
%   obj (CMSMSPepDeconv)
%       Updated instance
% Usage:
%   pairs = { ...
%       'runA.mgf', 'runB.mgf'; ...
%       'runC.mgf', 'runD.mgf' ...
%   };
%   strategy = PairwiseRunAlignStrategy(pairs);
%   align_options = struct('min_psm', 3, 'rt_sigma', 0.5, 'outlier_k', 3, ...
%       'outlier_method', 'mad', 'dead_time_min', 0.5);
%   obj = runImpAlignRequantLevel(obj, align_strategy, align_options);

if nargin < 2 || isempty(align_strategy)
    error('runImpAlignRequantLevel requires a CRunAlignStrategy instance.');
end
if nargin < 3
    align_options = struct();
end
if isempty(obj.m_filtered_res_file_path)
    error('FDR filtered result path is required for alignment.');
end

if ~isfolder(obj.m_outputDir)
    mkdir(obj.m_outputDir);
end

msms_res_path = COptionUtils.get(align_options, 'msms_res_path', obj.m_msms_res_path);
if isempty(msms_res_path)
    msms_res_path = fullfile(obj.m_outputDir, 'report_msms.txt');
end

% Read MSMS results
msms_result = CMS2ResultIO.read(msms_res_path);

% Initialize dataset IOs lazily
[obj, ~] = obj.ensureMs12DatasetIO();
[obj, mgf_created_here] = obj.ensureMgfDatasetIO();
if mgf_created_here
    cleanup_mgf = onCleanup(@() obj.m_cMgfDatasetIO.CloseAllFile());
end

% Initialize shared services lazily
[obj, ~] = obj.ensurePepProtService();
[obj, ~] = obj.ensureMsFileMapper();

% Build aligned RT range map
anchor_selector = CXICAlignAnchorSelector();
aligner = CXICAligner(anchor_selector, align_options);
align_cfg = obj.buildAlignmentPipelineConfig(aligner, align_strategy, align_options);
align_executor = CIMPXICAlignRequantExecutor(align_cfg);
pipeline = CPeptideLevelPipeline([], align_executor);

[pep_rtrange_map, align_report] = pipeline.buildAlignedRtRangeMap(...
    msms_result, obj.m_filtered_res_file_path, ...
    @(spectrum_list) obj.buildRawIdentManagerFromSpectrumList(spectrum_list));

alignment_report_path = COptionUtils.get(align_options, 'alignment_report_path', ...
    fullfile(obj.m_outputDir, 'report_alignment.txt'));
pipeline.writeAlignmentReport(align_report, alignment_report_path);

% Re-quantify using aligned RT ranges
output_path = COptionUtils.get(align_options, 'requant_output_path', ...
    fullfile(obj.m_outputDir, 'report_peptide_all_requant_aligned.txt'));
report = CIMPQuantReport();
print_progress = CPrintProgress(length(msms_result.Peptides));
proc_cfg = obj.buildIMPProcessingPipelineConfig();
proc_executor = CIMPProcessingExecutor(proc_cfg);
proc_pipeline = CPeptideLevelPipeline(proc_executor);

fprintf('Re-quantifying at peptide level (aligned)...')
for idx_psf = 1:length(msms_result.Peptides)
    print_progress = print_progress.update_show(idx_psf);
    peptide_sequence = msms_result.Peptides(idx_psf).peptide_sequence;
    cell_prot_name_pos = obj.CPepProtService.get_protein_name_pos(peptide_sequence);
    rawIdentManager = obj.buildRawIdentManagerFromSpectrumList( ...
        msms_result.Peptides(idx_psf).spectrum_list);
    block = proc_pipeline.requantifyPeptideBlock(cell_prot_name_pos, ...
        rawIdentManager, pep_rtrange_map);
    report = report.append_block(block);
end
print_progress.last_update();
fprintf('done.\n');

CIMPQuantResultIO.write(report, output_path);
end
