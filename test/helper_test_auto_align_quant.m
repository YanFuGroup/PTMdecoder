function helper_test_auto_align_quant(projectRootDir, testDataDir, outputDir)
% HELPER_TEST_AUTO_ALIGN_QUANT Runs the auto align quant test case
% Input:
%   projectRootDir (1 x 1 char/string) - root of the codebase
%   testDataDir (1 x 1 char/string) - directory containing test data (e.g., .../code/test/data)
%   outputDir (1 x 1 char/string) - directory to write outputs
% Output:
%   (none)

    % --- Configuration ---
    result_dir = fullfile(outputDir, 'auto_align_quant');
    input_dir = fullfile(testDataDir, 'auto_align_quant');
    
    ms1_tolerance.value = 10;
    ms1_tolerance.isppm = 1;
    resFilterThres = 0.1;
    alpha = 0.01;

    mod_file = fullfile(projectRootDir, 'modify.ini');
    fixed_mod = '';
    variable_mod = 'Acetyl[K];Methyl[K];Dimethyl[K];Trimethyl[K]';
    enzyme.name = 'trypsin';
    enzyme.limits = 14.015650;

    align_param_map = containers.Map();
    align_param_map('mod_file_path') = mod_file;
    align_param_map('fixed_mod') = fixed_mod;
    align_param_map('variable_mod') = variable_mod;
    align_param_map('spec_dir_path') = input_dir;
    align_param_map('ms1_tolerance') = ms1_tolerance;
    align_param_map('alpha') = alpha;
    align_param_map('result_filter_threshold') = resFilterThres;
    align_param_map('fasta_file_path') = fullfile(input_dir,'uniprotkb_human_histone_E_coli_comb_rever_czy_20231015.fasta');
    align_param_map('regular_express') = '>([^ ,]*)';
    align_param_map('output_dir_path') = result_dir;
    align_param_map('msms_res_path') = fullfile(input_dir, 'report_msms.txt');
    align_param_map('filtered_res_file_path') = fullfile(input_dir, 'filtered_result_mascot.txt');
    align_param_map('min_MSMS_num') = 1;

    align_cfg = CPeptideAlignRequantServiceConfig.fromParamMap(align_param_map);
    align_requant_service = CPeptideAlignRequantService(align_cfg);
    
    pairs = { ...
        'MCF7_JIB04_2_HISTONE_0723_HCDFT.mgf', 'MCF7_DMSO_2_HISTONE_0723_HCDFT.mgf'; ...
    };

    strategy = PairwiseRunAlignStrategy(pairs);
    align_options = struct('min_psm', 1, 'rt_sigma', 0.5, 'outlier_k', 3, ...
        'outlier_method', 'mad', 'dead_time_min', 0.5);
    align_requant_service.run(strategy, align_options);
end