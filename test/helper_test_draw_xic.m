function helper_test_draw_xic(projectRootDir, testDataDir, outputDir)
% HELPER_TEST_DRAW_XIC Runs the drawXIC test case
%
% Inputs:
%   projectRootDir - Root of the codebase
%   testDataDir    - Directory containing test data (e.g., .../code/test/data)
%   outputDir      - Directory to write outputs

    % --- Configuration ---
    
    % Data subdirectories
    draw_xic_data_dir = fullfile(testDataDir, 'draw_xic');
    msms_pep_site_dir = fullfile(testDataDir, 'msms_pep_site');
    
    % Input files
    pepSpec_path = fullfile(draw_xic_data_dir, 'pepSpecFile.txt');
    filtered_res_path = fullfile(msms_pep_site_dir, 'filtered_result_mascot.txt');
    checked_res_path = fullfile(draw_xic_data_dir, 'report_peptide_all.txt');
    msms_res_path = fullfile(draw_xic_data_dir, 'report_msms.txt');
    
    % Resource files
    modFile = fullfile(projectRootDir, 'modify.ini');
    fastaFile = fullfile(msms_pep_site_dir, 'uniprotkb_human_histone_E_coli_comb_rever_czy_20231015.fasta');
    
    parse_reg_exp = '>([^ ,]*)';

    % Parameters
    fixedMod = '';
    variableMod = 'Acetyl[K];Methyl[K];Dimethyl[K];Trimethyl[K]';

    enzyme.name = 'trypsin';
    enzyme.limits = 14.015650;
    
    lambda = 0.5;
    ms1_tolerance.value = 10;
    ms1_tolerance.isppm = 1;
    ms2_tolerance.value = 0.02;
    alpha = 0.01;
    resFilterThres = 0.1;
    
    ionTypes = [1,2];
    model = 1;
    method = 2;

    % --- Execution ---
    
    fprintf('Running test_draw_xic...\n');
    
    % Create processor object
    msms_pep_process = CMSMSPepDeconv(modFile, fixedMod, variableMod, ...
        msms_pep_site_dir, ms1_tolerance, ms2_tolerance.value, alpha, ...
        fastaFile, parse_reg_exp, pepSpec_path, model, method, lambda, ...
        resFilterThres, enzyme, outputDir, ionTypes, checked_res_path, ...
        msms_res_path, filtered_res_path);
    
    % Prepare maps
    [color_map, legend_map] = get_colormap_legendmap();
    
    % Run drawXIC
    msms_pep_process.drawXIC(outputDir, color_map, legend_map);
end

function [color_map, legend_map] = get_colormap_legendmap()
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
