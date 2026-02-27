function helper_test_draw_xic(projectRootDir, testDataDir, outputDir)
% HELPER_TEST_DRAW_XIC Runs the drawXIC test case
% Input:
%   projectRootDir (1 x 1 char/string) - root of the codebase
%   testDataDir (1 x 1 char/string) - directory containing test data (e.g., .../code/test/data)
%   outputDir (1 x 1 char/string) - directory to write outputs
% Output:
%   (none)

    % --- Configuration ---
    
    % Data subdirectories
    draw_xic_data_dir = fullfile(testDataDir, 'draw_xic');
    msms_pep_site_dir = fullfile(testDataDir, 'msms_pep_site');
    
    % Input files
    checked_res_path = fullfile(draw_xic_data_dir, 'report_peptide_all.txt');
    msms_res_path = fullfile(draw_xic_data_dir, 'report_msms.txt');
    
    % Resource files
    modFile = fullfile(projectRootDir, 'modify.ini');

    % Parameters
    fixedMod = '';
    variableMod = 'Acetyl[K];Methyl[K];Dimethyl[K];Trimethyl[K]';
    
    ms1_tolerance_value = 10;
    ms1_tolerance_type = 'PPM';
    alpha = 0.01;
    resFilterThres = 0.1;

    % --- Execution ---
    
    fprintf('Running test_draw_xic...\n');
    
    % Create service object
    draw_param_map = containers.Map();
    draw_param_map(CPTMdecoderWorkflowParamKeys.PARAM_MOD_FILE_PATH) = modFile;
    draw_param_map(CPTMdecoderWorkflowParamKeys.PARAM_FIXED_MOD) = fixedMod;
    draw_param_map(CPTMdecoderWorkflowParamKeys.PARAM_VARIABLE_MOD) = variableMod;
    draw_param_map(CPTMdecoderWorkflowParamKeys.PARAM_SPEC_DIR_PATH) = msms_pep_site_dir;
    draw_param_map(CPTMdecoderWorkflowParamKeys.PARAM_MS1_TOLERANCE_VALUE) = ms1_tolerance_value;
    draw_param_map(CPTMdecoderWorkflowParamKeys.PARAM_MS1_TOLERANCE_TYPE) = ms1_tolerance_type;
    draw_param_map(CPTMdecoderWorkflowParamKeys.PARAM_ALPHA) = alpha;
    draw_param_map(CPTMdecoderWorkflowParamKeys.PARAM_RESULT_FILTER_THRESHOLD) = resFilterThres;
    draw_param_map(CPTMdecoderWorkflowParamKeys.PARAM_OUTPUT_DIR_PATH) = outputDir;
    draw_param_map(CPTMdecoderWorkflowParamKeys.PARAM_CHECKED_PEPTIDES_RES_PATH) = checked_res_path;
    draw_param_map(CPTMdecoderWorkflowParamKeys.PARAM_MSMS_RES_PATH) = msms_res_path;
    draw_param_map(CPTMdecoderWorkflowParamKeys.PARAM_MIN_MSMS_NUM) = 1;

    draw_cfg = CXICDrawServiceConfig.fromParamMap(draw_param_map);
    draw_service = CXICDrawService(draw_cfg);
    
    % Prepare maps
    [color_map, legend_map] = get_colormap_legendmap();
    
    % Run drawXIC
    draw_service.run(outputDir, color_map, legend_map);
end

function [color_map, legend_map] = get_colormap_legendmap()
% GET_COLORMAP_LEGENDMAP Build color/legend maps for drawXIC
% Input:
%   (none)
% Output:
%   color_map (containers.Map) - map from peptide string to RGB color
%   legend_map (containers.Map) - map from peptide string to legend label
    % Target peptides
    peptides = {
        '_K{Trimethyl}SAPATGGVK{Dimethyl}KPHR_';
        '_K{Trimethyl}SAPATGGVKK{Dimethyl}PHR_'
    };

    % Define the color_map
    color_map = containers.Map();
    colors = [
        [0, 0.4470, 0.7410];
        [0.8500, 0.3250, 0.0980]
    ];

    % Assign the peptides and colors to the color_map
    for i = 1:length(peptides)
        color_map(peptides{i}) = colors(i, :);
    end

    % Define the legend_map
    legend_map = containers.Map();
    legend_strings = {
        'IMP1';
        'IMP2'
    };

    % Assign the peptides and legends to the legend_map
    for i = 1:length(peptides)
        legend_map(peptides{i}) = legend_strings{i};
    end
end
