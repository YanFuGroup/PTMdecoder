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
    
    ms1_tolerance_value = 10;
    ms1_tolerance_type = 'PPM';
    resFilterThres = 0.1;
    alpha = 0.01;

    mod_file = fullfile(projectRootDir, 'modify.ini');
    fixed_mod = '';
    variable_mod = 'Acetyl[K];Methyl[K];Dimethyl[K];Trimethyl[K]';

    align_param_map = containers.Map();
    align_param_map(CPTMdecoderWorkflowParamKeys.PARAM_MOD_FILE_PATH) = mod_file;
    align_param_map(CPTMdecoderWorkflowParamKeys.PARAM_FIXED_MOD) = fixed_mod;
    align_param_map(CPTMdecoderWorkflowParamKeys.PARAM_VARIABLE_MOD) = variable_mod;
    align_param_map(CPTMdecoderWorkflowParamKeys.PARAM_SPEC_DIR_PATH) = input_dir;
    align_param_map(CPTMdecoderWorkflowParamKeys.PARAM_MS1_TOLERANCE_VALUE) = ms1_tolerance_value;
    align_param_map(CPTMdecoderWorkflowParamKeys.PARAM_MS1_TOLERANCE_TYPE) = ms1_tolerance_type;
    align_param_map(CPTMdecoderWorkflowParamKeys.PARAM_ALPHA) = alpha;
    align_param_map(CPTMdecoderWorkflowParamKeys.PARAM_RESULT_FILTER_THRESHOLD) = resFilterThres;
    align_param_map(CPTMdecoderWorkflowParamKeys.PARAM_FASTA_FILE_PATH) = fullfile(input_dir,'uniprotkb_human_histone_E_coli_comb_rever_czy_20231015.fasta');
    align_param_map(CPTMdecoderWorkflowParamKeys.PARAM_REGULAR_EXPRESS) = '>([^ ,]*)';
    align_param_map(CPTMdecoderWorkflowParamKeys.PARAM_OUTPUT_DIR_PATH) = result_dir;
    align_param_map(CPTMdecoderWorkflowParamKeys.PARAM_MSMS_RES_PATH) = fullfile(input_dir, 'report_msms.txt');
    align_param_map(CPTMdecoderWorkflowParamKeys.PARAM_FILTERED_RES_FILE_PATH) = fullfile(input_dir, 'filtered_result_mascot.txt');
    align_param_map(CPTMdecoderWorkflowParamKeys.PARAM_MIN_MSMS_NUM) = 1;
    align_param_map(CPTMdecoderWorkflowParamKeys.PARAM_ALIGN_STRATEGY_TYPE) = 'pairwise';
    align_param_map(CPTMdecoderWorkflowParamKeys.PARAM_ALIGN_PAIR_NUM) = 1;
    align_param_map([CPTMdecoderWorkflowParamKeys.PARAM_PREFIX_ALIGN_PAIR, '1']) = ...
        'MCF7_JIB04_2_HISTONE_0723_HCDFT.mgf|MCF7_DMSO_2_HISTONE_0723_HCDFT.mgf';
    align_param_map(CPTMdecoderWorkflowParamKeys.PARAM_ALIGN_MIN_PSM) = 1;
    align_param_map(CPTMdecoderWorkflowParamKeys.PARAM_ALIGN_RT_SIGMA) = 0.5;
    align_param_map(CPTMdecoderWorkflowParamKeys.PARAM_ALIGN_OUTLIER_K) = 3;
    align_param_map(CPTMdecoderWorkflowParamKeys.PARAM_ALIGN_OUTLIER_METHOD) = 'mad';
    align_param_map(CPTMdecoderWorkflowParamKeys.PARAM_ALIGN_DEAD_TIME_MIN) = 0.5;

    align_cfg = CPeptideAlignRequantServiceConfig.fromParamMap(align_param_map);
    align_requant_service = CPeptideAlignRequantService(align_cfg);
    align_requant_service.run();
end