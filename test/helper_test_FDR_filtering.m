function helper_test_FDR_filtering(projectRootDir, testDataDir, outputDir)
% HELPER_TEST_FDR_FILTERING Runs the FDR filtering test case
% Input:
%   projectRootDir (1 x 1 char/string) - root of the codebase
%   testDataDir (1 x 1 char/string) - directory containing test data (e.g., .../code/test/data)
%   outputDir (1 x 1 char/string) - directory to write outputs
% Output:
%   (none)

    % --- Configuration ---
    result_dir = fullfile(outputDir, 'FDR_filtering');
    input_dir = fullfile(testDataDir, 'FDR_filtering');
    FDR_control_code_dir = fullfile(projectRootDir, 'FDR_control_generate_pep_spec_list');

    if ~exist(result_dir, 'dir')
        mkdir(result_dir);
    end

    addpath(FDR_control_code_dir);
    
    result = ReadDatResultFolder(input_dir);
    
    % Note: This needs to be modified to filter out results that do not conform to the experimental background.
    % For example, too many types of modifications, or modifications that should not occur at the C-terminus.
    result = FilterWithChemPrior(result);

    %% Judge Group
    DecoyTag = 'REVERSE_';
    GroupTag = {'_HUMAN'};

    [DecoyType,GroupType,~,scores,numrst,I] = JudgeGroup(result,'Protein',DecoyTag,GroupTag);
    %% Compute FDR
    fdrthres = 0.01;
    [~,Iid,~,~] = ComputeFDR(DecoyType,GroupType,scores,numrst,I,fdrthres);
    fprintf('The size of the PSM in target group is %d.\n',sum(~GroupType));
    fprintf('The size of filtered PSM is %d.\n',length(Iid.SF));

    %% select using some criteria and write to file
    % Write all grouped result
    group_filtered_result = result(I(~GroupType));
    grouped_filename = fullfile(result_dir,'group_result_mascot.txt');
    write_mascot_result_table(group_filtered_result, grouped_filename);

    % Write all filtered result
    FDR_filtered_result = result(Iid.SF);
    filtered_fileName = fullfile(result_dir,'filtered_result_mascot.txt');
    write_mascot_result_table(FDR_filtered_result, filtered_fileName);

    % select the unique identified spectra
    mgf_scan_strings = strcat({FDR_filtered_result.DatasetName}',' ',{FDR_filtered_result.Scan}');
    is_result_selected = false(length(FDR_filtered_result),1);
    for i_res = 1:length(FDR_filtered_result)
        if sum(strcmp(mgf_scan_strings(i_res),mgf_scan_strings)) == 1
            % set the marker to retain this unique element
            is_result_selected(i_res) = true;
        end
    end
    unique_filtered_result = FDR_filtered_result(is_result_selected);

    % select the modified peptides without Oxidation
    is_result_selected = false(length(unique_filtered_result),1);
    for i_final = 1:length(unique_filtered_result)
        if ~isequal(unique_filtered_result(i_final).modification,'-')...
                && ~contains(unique_filtered_result(i_final).modification,'Oxidation')
            is_result_selected(i_final) = true;
        end
    end
    final_result = unique_filtered_result(is_result_selected);

    % sort the results by peptide sequence
    [~,sorted_idx] = sort({final_result.peptide});
    output_res = final_result(sorted_idx);

    % write pep_spec to files
    fileName = fullfile(result_dir,'pepSpecFile.txt');
    write_peptide_spectra_list_file(output_res, fileName);

    rmpath(FDR_control_code_dir);
end