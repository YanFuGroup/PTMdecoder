function obj = quant_norm_pep_from_fdr(obj, peptide_list, prot_list, filtered_res_file_path, output_file_name)
% Quantify normalization peptides from FDR results (single experiment)
% Input:
%   obj (CMSMSPepDeconv)
%       Processor instance
%   peptide_list (K x 1 cellstr/string)
%       Target peptide sequences
%   prot_list (K x 1 cellstr/string)
%       Corresponding protein names
%   filtered_res_file_path (1 x 1 char/string)
%       Path to filtered PSM result file
%   output_file_name (1 x 1 char/string)
%       Output report file name

if nargin < 5 || isempty(output_file_name)
    output_file_name = 'peptide4normalization.txt';
end

if length(peptide_list) ~= length(prot_list)
    error('peptide_list and prot_list must have the same length');
end

if ~exist(obj.m_outputDir, 'dir')
    mkdir(obj.m_outputDir);
end

% Indexing the ms1 and ms2
ms12DatasetIO = CMS12DatasetIO(obj.m_specPath, obj.m_ms1_tolerance);
ms12DatasetIO.SetMap();

% Initialize raw managers
pep_quant = cell(length(peptide_list), 1);
for i_list = 1:length(peptide_list)
    pep_quant{i_list} = CIMPRawIdentManager();
end

% Process the filtered result file
fprintf('Reading %s...', obj.m_specPath);
pep_quant = obj.readFdrPeptides(filtered_res_file_path, ms12DatasetIO, peptide_list, pep_quant);
fprintf('done.\n');

% Run quantification
pipeline = CIMPProcessingPipeline(ms12DatasetIO, obj.m_ms1_tolerance, 1, obj.m_alpha, obj.m_resFilterThres);
stats_cleanup = onCleanup(@() CIMPQuantifier.rt_sorted_stats('flush', fullfile(obj.m_outputDir, 'rt_sorted_stats.mat')));
report = CIMPQuantReport();
fprintf('Quantifying %s...', obj.m_specPath);
for i_list = 1:length(peptide_list)
    block = pipeline.quantifyPeptideBlock({prot_list{i_list}, -1}, pep_quant{i_list});
    report = report.append_block(block);
end
fprintf('done.\n');

CIMPQuantResultIO.write(report, fullfile(obj.m_outputDir, output_file_name));
end
