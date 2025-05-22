clear
%% read Mascot(*.Dat) path
res_path = 'D:\research\project\Mixspec_code\paper_data\MCF7_JIB04\mascot_dat';
work_dir_res = 'D:\research\project\Mixspec_code\paper_data\MCF7_JIB04\result';
work_dir_res_top1 = 'D:\research\project\Mixspec_code\paper_data\MCF7_JIB04\result_top1';
experimentNames = {'DMSO_1','DMSO_2','JIB04_1','JIB04_2'};
for idx = 1:length(experimentNames)
    result = ReadDatResultFolder(fullfile(res_path,experimentNames{idx}));
    
    % Note: This needs to be modified to filter out results that do not conform to the experimental background.
    % For example, too many types of modifications, or modifications that should not occur at the C-terminus.
    result = FilterWithChemPrior(result);

    %% Judge Group
    DecoyTag = 'REVERSE_';
    GroupTag = {'_HUMAN'};

    [DecoyType,GroupType,OriginalPeptide,scores,numrst,I] = JudgeGroup(result,'Protein',DecoyTag,GroupTag);
    %% Compute FDR
    fdrthres = 0.01;
    [FDR,Iid,threshold,finalFDR] = ComputeFDR(DecoyType,GroupType,scores,numrst,I,fdrthres);
    fprintf('%s: The size of the PSM in target group is %d.\n',experimentNames{idx},sum(~GroupType));
    fprintf('%s: The size of filtered PSM is %d.\n',experimentNames{idx},length(Iid.SF));

    %% select using some criteria and write to file
    % Write all grouped result
    group_filtered_result = result(I(~GroupType));
    grouped_filename = fullfile(work_dir_res,experimentNames{idx},'group_result_mascot.txt');
    write_mascot_result_table(group_filtered_result, grouped_filename);
    grouped_filename = fullfile(work_dir_res_top1,experimentNames{idx},'group_result_mascot.txt');
    write_mascot_result_table(group_filtered_result, grouped_filename);

    % Write all filtered result
    FDR_filtered_result = result(Iid.SF);
    filtered_fileName = fullfile(work_dir_res,experimentNames{idx},'filtered_result_mascot.txt');
    write_mascot_result_table(FDR_filtered_result, filtered_fileName);
    filtered_fileName = fullfile(work_dir_res_top1,experimentNames{idx},'filtered_result_mascot.txt');
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
    fileName = fullfile(work_dir_res,experimentNames{idx},'pepSpecFile.txt');
    write_peptide_spectra_list_file(output_res, fileName);
    % write the top scoring results to file in report_msms.txt format
    fileName = fullfile(work_dir_res_top1,experimentNames{idx},'report_msms_top1.txt');
    write_report_msms_top1(output_res, fileName);

end