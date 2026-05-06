function tests = test_CTraceabilityReportService
% TEST_CTRACEABILITYREPORTSERVICE Verify traceability sidecar report generation.

tests = functiontests(localfunctions);
end


function testTraceReportsUseKeysAsFirstColumnsAndMgfCharge(testCase)
% Verify both trace files and charge lookup by exact MGF TITLE.
workDir = tempname;
mkdir(workDir);

specDir = fullfile(workDir, 'spectra');
mkdir(specDir);
writeLines(fullfile(specDir, 'sample.mgf'), {
    'BEGIN IONS'
    'TITLE=RawFile.1234.1234.2.0.dta'
    'PEPMASS=500.25'
    'CHARGE=2+'
    '100.1 10'
    'END IONS'
});

msmsPath = fullfile(workDir, 'report_msms.txt');
writeLines(msmsPath, {
    'P	AK'
    'S	sample.mgf	RawFile.1234.1234.2.0.dta'
    'AK{Acetyl}	42.5	support=0.75	mad=0.125'
});

peptidePath = fullfile(workDir, 'report_peptide_all_requant_aligned.txt');
writeLines(peptidePath, {
    'Protein_name,Peptide_start_position_on_protein;'
    '*	IMP	Charge	Dataset	Mass_center	Low_mass_bound	High_mass_bound	Peak_area'
    '@	RT_start	RT_end	Proportion	Check_label'
    'P1,10;P2,30;'
    '*	_AK{Acetyl}_	+2	sample.mgf	500.25	499.95	500.55	100'
});

cfg = makeTraceConfig(specDir, msmsPath, peptidePath, workDir);
service = CTraceabilityReportService(cfg);
cleanupAll = onCleanup(@() cleanupServiceAndDir(service, workDir));
service.run();

peptideMsmsLines = readLines(cfg.output_trace_peptide_msms_path);
peptideMsmsHeader = strsplit(peptideMsmsLines{1}, sprintf('\t'));
peptideMsmsValues = strsplit(peptideMsmsLines{2}, sprintf('\t'));
testCase.verifyEqual(peptideMsmsHeader{1}, 'imp_msms_key');
testCase.verifyEqual(peptideMsmsValues{1}, 'sample.mgf|AK{Acetyl}|+2');
testCase.verifyEqual(peptideMsmsValues{2}, 'sample.mgf|RawFile.1234.1234.2.0.dta');
testCase.verifyEqual(peptideMsmsValues{5}, '2');
testCase.verifyEqual(peptideMsmsValues{8}, '42.5');
testCase.verifyEqual(peptideMsmsValues{9}, '0.75');
testCase.verifyEqual(peptideMsmsValues{10}, '0.125');

sitePeptideLines = readLines(cfg.output_trace_site_peptide_path);
sitePeptideHeader = strsplit(sitePeptideLines{1}, sprintf('\t'));
sitePeptideValues = strsplit(sitePeptideLines{2}, sprintf('\t'));
testCase.verifyEqual(sitePeptideHeader{1}, 'site_peptide_key');
testCase.verifyEqual(sitePeptideHeader{2}, 'site_name');
testCase.verifyEqual(sitePeptideHeader{8}, 'protein_site_positions');
testCase.verifyEqual(sitePeptideValues{1}, 'P1,11;P2,31; K_ac|sample.mgf|_AK{Acetyl}_|+2');
testCase.verifyEqual(sitePeptideValues{2}, 'P1,11;P2,31; K_ac');
testCase.verifyEqual(sitePeptideValues{5}, '2');
testCase.verifyEqual(sitePeptideValues{6}, 'AK');
testCase.verifyEqual(sitePeptideValues{7}, 'P1,10;P2,30;');
testCase.verifyEqual(sitePeptideValues{8}, 'P1,11;P2,31;');
testCase.verifyEqual(sitePeptideValues{15}, '100');
end


function testMissingSpectrumChargeThrowsWithoutFallback(testCase)
% Verify spectrum_name must directly match the charge map key.
workDir = tempname;
mkdir(workDir);

specDir = fullfile(workDir, 'spectra');
mkdir(specDir);
writeLines(fullfile(specDir, 'sample.mgf'), {
    'BEGIN IONS'
    'TITLE=UnparseableTitle'
    'PEPMASS=500.25'
    'CHARGE=2+'
    'END IONS'
});

msmsPath = fullfile(workDir, 'report_msms.txt');
writeLines(msmsPath, {
    'P	AK'
    'S	sample.mgf	1234'
    'AK{Acetyl}	42.5'
});

peptidePath = fullfile(workDir, 'report_peptide_all_requant_aligned.txt');
writeLines(peptidePath, {
    'Protein_name,Peptide_start_position_on_protein;'
    '*	IMP	Charge	Dataset	Mass_center	Low_mass_bound	High_mass_bound	Peak_area'
    '@	RT_start	RT_end	Proportion	Check_label'
    'P1,10;'
    '*	_AK{Acetyl}_	+2	sample.mgf	500.25	499.95	500.55	100'
});

cfg = makeTraceConfig(specDir, msmsPath, peptidePath, workDir);
service = CTraceabilityReportService(cfg);
cleanupAll = onCleanup(@() cleanupServiceAndDir(service, workDir));
testCase.verifyError(@() service.run(), 'CTraceabilityReportService:MissingSpectrumCharge');
end


function testSitePeptideTraceKeepsChargesAndAggregates(testCase)
% Verify site records keep charge as an attribute and can reproduce site sums.
workDir = tempname;
mkdir(workDir);

specDir = fullfile(workDir, 'spectra');
mkdir(specDir);
writeLines(fullfile(specDir, 'sample.mgf'), {
    'BEGIN IONS'
    'TITLE=SpecA'
    'PEPMASS=500.25'
    'CHARGE=2+'
    'END IONS'
});

msmsPath = fullfile(workDir, 'report_msms.txt');
writeLines(msmsPath, {
    'P	AK'
    'S	sample.mgf	SpecA'
    'AK{Acetyl}	1'
});

peptidePath = fullfile(workDir, 'report_peptide_all_requant_aligned.txt');
writeLines(peptidePath, {
    'Protein_name,Peptide_start_position_on_protein;'
    '*	IMP	Charge	Dataset	Mass_center	Low_mass_bound	High_mass_bound	Peak_area'
    '@	RT_start	RT_end	Proportion	Check_label'
    'P1,10;'
    '*	_AK{Acetyl}_	+2	sample.mgf	500.25	499.95	500.55	100'
    '*	_AK{Acetyl}_	+3	sample.mgf	333.25	332.95	333.55	50'
});

cfg = makeTraceConfig(specDir, msmsPath, peptidePath, workDir);
service = CTraceabilityReportService(cfg);
cleanupAll = onCleanup(@() cleanupServiceAndDir(service, workDir));
service.run();

sitePeptideLines = readLines(cfg.output_trace_site_peptide_path);
testCase.verifyEqual(numel(sitePeptideLines), 3);

firstValues = strsplit(sitePeptideLines{2}, sprintf('\t'));
secondValues = strsplit(sitePeptideLines{3}, sprintf('\t'));
testCase.verifyEqual(firstValues{1}, 'P1,11; K_ac|sample.mgf|_AK{Acetyl}_|+2');
testCase.verifyEqual(secondValues{1}, 'P1,11; K_ac|sample.mgf|_AK{Acetyl}_|+3');

areaSum = str2double(firstValues{15}) + str2double(secondValues{15});
testCase.verifyEqual(areaSum, 150);
end


function testSitePeptideTraceUsesConfiguredColumnIdxs(testCase)
% Verify peptide-level parsing uses cfg.column_idxs instead of fixed columns.
workDir = tempname;
mkdir(workDir);

specDir = fullfile(workDir, 'spectra');
mkdir(specDir);
writeLines(fullfile(specDir, 'sample.mgf'), {
    'BEGIN IONS'
    'TITLE=SpecA'
    'PEPMASS=500.25'
    'CHARGE=2+'
    'END IONS'
});

msmsPath = fullfile(workDir, 'report_msms.txt');
writeLines(msmsPath, {
    'P	AK'
    'S	sample.mgf	SpecA'
    'AK{Acetyl}	1'
});

peptidePath = fullfile(workDir, 'report_peptide_custom_columns.txt');
writeLines(peptidePath, {
    'Protein_name,Peptide_start_position_on_protein;'
    '*	Dataset	IMP	Peak_area	Charge	Mass_center	Low_mass_bound	High_mass_bound'
    '@	RT_start	RT_end	Proportion	Check_label'
    'P1,10;'
    '*	sample.mgf	_AK{Acetyl}_	100	+2	500.25	499.95	500.55'
});

cfg = makeTraceConfig(specDir, msmsPath, peptidePath, workDir);
cfg.column_idxs = struct('icol_seq', 3, 'icol_charge', 5, 'icol_dataset', 2, ...
    'icol_mass_center', 6, 'icol_low_mz_bound', 7, 'icol_high_mz_bound', 8, 'icol_auc', 4);

service = CTraceabilityReportService(cfg);
cleanupAll = onCleanup(@() cleanupServiceAndDir(service, workDir));
service.run();

sitePeptideLines = readLines(cfg.output_trace_site_peptide_path);
sitePeptideValues = strsplit(sitePeptideLines{2}, sprintf('\t'));
testCase.verifyEqual(sitePeptideValues{1}, 'P1,11; K_ac|sample.mgf|_AK{Acetyl}_|+2');
testCase.verifyEqual(sitePeptideValues{12}, '500.25');
testCase.verifyEqual(sitePeptideValues{15}, '100');
end


function testSitePeptideTraceCanExcludeInitialM(testCase)
% Verify trace site positions honor site_position_count_initial_m=false.
workDir = tempname;
mkdir(workDir);

specDir = fullfile(workDir, 'spectra');
mkdir(specDir);
writeLines(fullfile(specDir, 'sample.mgf'), {
    'BEGIN IONS'
    'TITLE=SpecA'
    'PEPMASS=500.25'
    'CHARGE=2+'
    'END IONS'
});

msmsPath = fullfile(workDir, 'report_msms.txt');
writeLines(msmsPath, {
    'P	AK'
    'S	sample.mgf	SpecA'
    'AK{Acetyl}	1'
});

peptidePath = fullfile(workDir, 'report_peptide_all_requant_aligned.txt');
writeLines(peptidePath, {
    'Protein_name,Peptide_start_position_on_protein;'
    '*	IMP	Charge	Dataset	Mass_center	Low_mass_bound	High_mass_bound	Peak_area'
    '@	RT_start	RT_end	Proportion	Check_label'
    'P1,10;P2,30;'
    '*	_AK{Acetyl}_	+2	sample.mgf	500.25	499.95	500.55	100'
});

cfg = makeTraceConfig(specDir, msmsPath, peptidePath, workDir);
cfg.site_position_count_initial_m = false;
service = CTraceabilityReportService(cfg);
cleanupAll = onCleanup(@() cleanupServiceAndDir(service, workDir));
service.run();

sitePeptideLines = readLines(cfg.output_trace_site_peptide_path);
sitePeptideValues = strsplit(sitePeptideLines{2}, sprintf('\t'));
testCase.verifyEqual(sitePeptideValues{1}, 'P1,10;P2,30; K_ac|sample.mgf|_AK{Acetyl}_|+2');
testCase.verifyEqual(sitePeptideValues{2}, 'P1,10;P2,30; K_ac');
testCase.verifyEqual(sitePeptideValues{7}, 'P1,10;P2,30;');
testCase.verifyEqual(sitePeptideValues{8}, 'P1,10;P2,30;');
end


function testSitePeptideTraceComputesTerminalPositions(testCase)
% Verify N-term and C-term site coordinates are traced on all proteins.
workDir = tempname;
mkdir(workDir);

specDir = fullfile(workDir, 'spectra');
mkdir(specDir);
writeLines(fullfile(specDir, 'sample.mgf'), {
    'BEGIN IONS'
    'TITLE=SpecA'
    'PEPMASS=500.25'
    'CHARGE=2+'
    'END IONS'
});

msmsPath = fullfile(workDir, 'report_msms.txt');
writeLines(msmsPath, {
    'P	AK'
    'S	sample.mgf	SpecA'
    'AK{Acetyl}	1'
});

peptidePath = fullfile(workDir, 'report_peptide_all_requant_aligned.txt');
writeLines(peptidePath, {
    'Protein_name,Peptide_start_position_on_protein;'
    '*	IMP	Charge	Dataset	Mass_center	Low_mass_bound	High_mass_bound	Peak_area'
    '@	RT_start	RT_end	Proportion	Check_label'
    'P1,10;P2,30;'
    '*	_{Acetyl}AK_	+2	sample.mgf	500.25	499.95	500.55	100'
    '*	_AK_{Acetyl}	+2	sample.mgf	500.25	499.95	500.55	50'
});

cfg = makeTraceConfig(specDir, msmsPath, peptidePath, workDir);
service = CTraceabilityReportService(cfg);
cleanupAll = onCleanup(@() cleanupServiceAndDir(service, workDir));
service.run();

sitePeptideLines = readLines(cfg.output_trace_site_peptide_path);
nTermValues = strsplit(sitePeptideLines{2}, sprintf('\t'));
cTermValues = strsplit(sitePeptideLines{3}, sprintf('\t'));
testCase.verifyEqual(nTermValues{2}, 'P1,9;P2,29; N-term_ac');
testCase.verifyEqual(nTermValues{8}, 'P1,9;P2,29;');
testCase.verifyEqual(nTermValues{9}, 'N-term');
testCase.verifyEqual(cTermValues{2}, 'P1,12;P2,32; C-term_ac');
testCase.verifyEqual(cTermValues{8}, 'P1,12;P2,32;');
testCase.verifyEqual(cTermValues{9}, 'C-term');
end


function testConfigDefaultsAndWorkflowStage(testCase)
% Verify default trace output paths and workflow stage registration.
workDir = tempname;
mkdir(workDir);
cleanupWork = onCleanup(@() cleanupTempDir(workDir));

paramPath = fullfile(workDir, 'params.txt');
writeLines(paramPath, {
    'traceability_report_on=1'
    ['spec_dir_path=', fullfile(workDir, 'spectra')]
    ['msms_res_path=', fullfile(workDir, 'report_msms.txt')]
    ['pep_level_file_path=', fullfile(workDir, 'report_peptide_all_requant_aligned.txt')]
    ['output_dir_path=', workDir]
    'mod_name_abbr_num=1'
    'mod_name_abbr_1=Acetyl>ac'
    'ignore_strings_site_level=""'
    'log_enabled=0'
});

workflowCfg = CPTMdecoderWorkflowConfig.fromParamFile(paramPath);
testCase.verifyEqual(numel(workflowCfg.stages), 1);
testCase.verifyEqual(workflowCfg.stages{1}.name, CPTMdecoderWorkflowConfig.STAGE_TRACEABILITY_REPORT);

traceCfg = workflowCfg.stages{1}.config;
testCase.verifyEqual(traceCfg.output_trace_peptide_msms_path, ...
    fullfile(workDir, 'report_trace_peptide_msms.txt'));
testCase.verifyEqual(traceCfg.output_trace_site_peptide_path, ...
    fullfile(workDir, 'report_trace_site_peptide.txt'));
testCase.verifyTrue(traceCfg.site_position_count_initial_m);
end


function testConfigCanDisableInitialMCounting(testCase)
% Verify traceability config can disable initial-M site numbering.
workDir = tempname;
mkdir(workDir);
cleanupWork = onCleanup(@() cleanupTempDir(workDir));

paramPath = fullfile(workDir, 'params.txt');
writeLines(paramPath, {
    'traceability_report_on=1'
    ['spec_dir_path=', fullfile(workDir, 'spectra')]
    ['msms_res_path=', fullfile(workDir, 'report_msms.txt')]
    ['pep_level_file_path=', fullfile(workDir, 'report_peptide_all_requant_aligned.txt')]
    ['output_dir_path=', workDir]
    'mod_name_abbr_num=1'
    'mod_name_abbr_1=Acetyl>ac'
    'ignore_strings_site_level=""'
    'site_position_count_initial_m=0'
    'log_enabled=0'
});

workflowCfg = CPTMdecoderWorkflowConfig.fromParamFile(paramPath);
traceCfg = workflowCfg.stages{1}.config;
testCase.verifyFalse(traceCfg.site_position_count_initial_m);
end


function cfg = makeTraceConfig(specDir, msmsPath, peptidePath, outputDir)
% Build a minimal traceability config for tests.
modMap = containers.Map('KeyType', 'char', 'ValueType', 'char');
modMap('Acetyl') = 'ac';

cfg = struct();
cfg.spec_dir_path = specDir;
cfg.msms_res_path = msmsPath;
cfg.pep_level_file_path = peptidePath;
cfg.output_trace_peptide_msms_path = fullfile(outputDir, 'report_trace_peptide_msms.txt');
cfg.output_trace_site_peptide_path = fullfile(outputDir, 'report_trace_site_peptide.txt');
cfg.mod_name_abbr = modMap;
cfg.ignore_strings = {};
cfg.site_position_count_initial_m = true;
cfg.column_idxs = struct('icol_seq', 2, 'icol_charge', 3, 'icol_dataset', 4, ...
    'icol_mass_center', 5, 'icol_low_mz_bound', 6, 'icol_high_mz_bound', 7, 'icol_auc', 8);
end


function writeLines(filePath, lines)
% Write text lines into a file with native newlines.
fid = fopen(filePath, 'w');
if fid < 0
    error('test_CTraceabilityReportService:FileOpenFailed', ...
        'Cannot open file: %s', filePath);
end
cleanupFile = onCleanup(@() fclose(fid));
for idx = 1:numel(lines)
    fprintf(fid, '%s\n', lines{idx});
end
end


function lines = readLines(filePath)
% Read nonempty lines from a text file.
content = fileread(filePath);
rawLines = regexp(content, '\r\n|\n|\r', 'split');
lines = rawLines(~cellfun('isempty', rawLines));
end


function cleanupServiceAndDir(service, workDir)
% Release service-owned files before deleting the temporary directory.
if ~isempty(service)
    delete(service);
end
cleanupTempDir(workDir);
end


function cleanupTempDir(workDir)
% Delete a temporary directory if it still exists.
if isfolder(workDir)
    rmdir(workDir, 's');
end
end
