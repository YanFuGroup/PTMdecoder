function obj = quant_IMP(obj)
% Run the quantification at peptide level
% Read the results of PSM level and quantify at peptide level
% Input:
%   obj (CMSMSPepDeconv)
%       Processor instance
% Output:
%   obj (CMSMSPepDeconv)
%       Updated instance

% Load protein sequences from fasta file
if isempty(obj.CPepProtService)
    obj.CPepProtService = CPepProtService(obj.m_fastaFile, obj.m_regular_express, obj.m_filtered_res_file_path);
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

if ~isfolder(obj.m_outputDir)
    mkdir(obj.m_outputDir);
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
report = CIMPQuantReport();

% Read and process
msms_reader = CMSMSResReader();
msms_result = msms_reader.read_from_msms_res_file(each_PSM_results_path);
print_progress = CPrintProgress(length(msms_result.Peptides));
pipeline = CIMPProcessingPipeline(obj.m_cMs12DatasetIO, obj.m_ms1_tolerance, obj.m_min_MSMS_num, obj.m_alpha, obj.m_resFilterThres);
stats_cleanup = onCleanup(@() CIMPQuantifier.rt_sorted_stats('flush', fullfile(obj.m_outputDir, 'rt_sorted_stats.mat')));

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

if need_release_mgf_index
    obj.m_cMgfDatasetIO.CloseAllFile();
end
end