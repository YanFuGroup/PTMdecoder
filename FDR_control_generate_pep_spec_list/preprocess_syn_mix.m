clear
%% read Mascot(*.Dat) path
res_path = 'D:\research\project\Mixspec_code\paper_data\syn_mix\mascot_dat';
work_dir_res = 'D:\research\project\Mixspec_code\paper_data\syn_mix\result';
experimentNames = {'mix1', 'mix2', 'mix3', 'mix4', 'mix5', 'mix6', 'mix7', 'mix8', ...
    'mix9', 'mix10', 'mix11', 'mix12', 'mix13', 'mix14', 'mix15', 'mix16', 'mix17', 'mix18'};
for idx = 1:length(experimentNames)
    result = ReadDatResultFolder(fullfile(res_path,experimentNames{idx}));
    % result = FilterWithChemPrior(result); % Do not need to filter with chemical prior

    %% Judge Group
    DecoyTag = 'REVERSE_';
    TagType = 'Protein';
    GroupTag = {'_HUMAN'};

    [DecoyType,GroupType,OriginalPeptide,scores,numrst,I] = JudgeGroup(result,TagType,DecoyTag,GroupTag);
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

    % Write all filtered result
    FDR_filtered_result = result(Iid.SF);
    filtered_fileName = fullfile(work_dir_res,experimentNames{idx},'filtered_result_mascot.txt');
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

    % select the target peptides GKGGKGLGKGGAKR
    is_result_selected = false(length(unique_filtered_result),1);
    for i_final = 1:length(unique_filtered_result)
        if isequal(unique_filtered_result(i_final).peptide,'GKGGKGLGKGGAKR') && ...
                abs(unique_filtered_result(i_final).Calc_neutral_pepmass-1535.878357)<0.0001
            is_result_selected(i_final) = true;
        end
    end
    final_result = unique_filtered_result(is_result_selected);

    % sort the results by peptide sequences
    [~,sorted_idx] = sort({final_result.peptide});
    output_res = final_result(sorted_idx);

    % write pep_spec to files
    fileName = fullfile(work_dir_res,experimentNames{idx},'pepSpecFile.txt');
    write_peptide_spectra_list_file(output_res, fileName);

end