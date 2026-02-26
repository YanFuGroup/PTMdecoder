function obj = runImpQuantLevel(obj)
% Run the quantification at peptide level
% Read the results of PSM level and quantify at peptide level
% Input:
%   obj (CMSMSPepDeconv)
%       Processor instance
% Output:
%   obj (CMSMSPepDeconv)
%       Updated instance

[obj, ~] = obj.ensurePepProtService();

obj = obj.ensureMgfDatasetIO();

[obj, ~] = obj.ensureMs12DatasetIO();

CPathResolver.ensureDir(obj.m_outputDir);

% Check the report_msms.txt file
each_PSM_results_path = CPathResolver.resolveFilePath(obj.m_outputDir, ...
	'report_msms.txt', obj.m_msms_res_path);

% Check and create a new output file
% TODO: Set the output path for peptide level results
each_peptide_results_path = CPathResolver.resolveFilePath(obj.m_outputDir, 'report_peptide_all.txt', '');
report = CIMPQuantReport();

% Read and process
msms_result = CMS2ResultIO.read(each_PSM_results_path);
print_progress = CPrintProgress(length(msms_result.Peptides));
pipeline_cfg = obj.buildIMPProcessingPipelineConfig();
executor = CIMPProcessingExecutor(pipeline_cfg);
pipeline = CPeptideLevelPipeline(executor);
stats_cleanup = onCleanup(@() CIMPQuantStats.rt_sorted_stats('flush', CPathResolver.resolveFilePath(obj.m_outputDir, 'rt_sorted_stats.mat', '')));

fprintf('Quantifying at peptide level...')
for idx_psf = 1:length(msms_result.Peptides)
	% Show progress
	print_progress = print_progress.update_show(idx_psf);
	% Get the peptide sequence
	peptide_sequence = msms_result.Peptides(idx_psf).peptide_sequence;
	% Get the protein name and position
	cell_prot_name_pos = obj.CPepProtService.get_protein_name_pos(peptide_sequence);
	rawIdentManager = obj.buildRawIdentManagerFromSpectrumList(msms_result.Peptides(idx_psf).spectrum_list);
	% Run gather
	block = pipeline.quantifyPeptideBlock(cell_prot_name_pos, rawIdentManager);
	report = report.append_block(block);
end
print_progress.last_update();
fprintf('done.\n');

CIMPQuantResultIO.write(report, each_peptide_results_path);
end
