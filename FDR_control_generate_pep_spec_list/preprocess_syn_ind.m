clear
%% read Mascot(*.Dat) path
res_path = 'D:\research\project\Mixspec_code\paper_data\syn_ind\mascot_dat';
workspacePath = 'D:\research\project\Mixspec_code\paper_data\syn_ind\result';
experimentNames = {'ind1','ind2','ind3','ind4','ind5','ind6','ind7','ind8','ind9','ind10', ...
    'ind11','ind12','ind13','ind14','ind15','ind16'};
for idx = 1:length(experimentNames)
    result = ReadDatResultFolder(fullfile(res_path,experimentNames{idx}));
    
    % Some prior chemical knowledge, such as there should be no more than 2
    %  kinds of modifications on a peptides.
%     result = FilterWithChemPrior(result);

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
    grouped_filename = fullfile(workspacePath,experimentNames{idx},'group_result_mascot.txt');
    write_mascot_result_table(group_filtered_result, grouped_filename);

    % Write all filtered result
    FDR_filtered_result = result(Iid.SF);
    filtered_fileName = fullfile(workspacePath,experimentNames{idx},'filtered_result_mascot.txt');
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
                isequal(unique_filtered_result(i_final).Charge, '2') && ...
                abs(unique_filtered_result(i_final).Calc_neutral_pepmass-1535.878357)<0.0001 && ...
                unique_filtered_result(i_final).num_match_ions >= 8
            is_result_selected(i_final) = true;
        end
    end
    final_result = unique_filtered_result(is_result_selected);

    % sort the results by peptide sequence
    [~,sorted_idx] = sort({final_result.peptide});
    output_res = final_result(sorted_idx);

    % write pep_spec to files
    fileName = fullfile(workspacePath,experimentNames{idx},'pepSpecFile.txt');
    write_peptide_spectra_list_file(output_res, fileName);
end