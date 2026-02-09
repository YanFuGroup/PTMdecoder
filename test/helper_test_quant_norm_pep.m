function helper_test_quant_norm_pep(projectRootDir, testDataDir, outputDir)
% HELPER_TEST_QUANT_NORM_PEP Runs the quant normalization peptide test case
% Input:
%   projectRootDir (1 x 1 char/string) - root of the codebase
%   testDataDir (1 x 1 char/string) - directory containing test data (e.g., .../code/test/data)
%   outputDir (1 x 1 char/string) - directory to write outputs
% Output:
%   (none)

    % --- Configuration ---
    result_dir = fullfile(outputDir, 'requant_norm_pep');
    input_dir = fullfile(testDataDir, 'requant_norm_pep');
    
    peptide_list = {
        'AGLQFPVGR','VTIAQGGVLPNIQAVLLPK',...
        'LLLPGELAK','EIQTAVR',...
        'STELLIR','EIAQDFK',...
        'TLYGFGG','VFLENVIR','ISGLIYEETR'
    };
    prot_list = {
        'H2A','H2A',...
        'H2B','H2B',...
        'H3','H3',...
        'H4','H4','H4'
    };
    experimentNames = {''};
    ms1_tolerance.value = 10;
    ms1_tolerance.isppm = 1;
    resFilterThres = 0.1;
    alpha = 0.01;

    normPepQuant = CPepNormalization(peptide_list, prot_list, result_dir, ...
    input_dir, experimentNames, ms1_tolerance, resFilterThres, alpha, fullfile(input_dir, 'filtered_result_mascot.txt'));
    normPepQuant.runQuant()
    normPepQuant.runRequant(fullfile(input_dir, 'peptide4normalization_checked.txt'));
end