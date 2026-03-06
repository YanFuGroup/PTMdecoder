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
    [~,~,group_idx] = unique(mgf_scan_strings,'stable');
    group_counts = accumarray(group_idx,1);
    is_result_selected = (group_counts(group_idx) == 1);
    unique_filtered_result = FDR_filtered_result(is_result_selected);

    % select the target peptides GKGGKGLGKGGAKR
    peptides = {unique_filtered_result.peptide}';
    charges = {unique_filtered_result.Charge}';
    calc_neutral_pepmass = [unique_filtered_result.Calc_neutral_pepmass]';
    num_match_ions = [unique_filtered_result.num_match_ions]';
    is_result_selected = strcmp(peptides,'GKGGKGLGKGGAKR') ...
        & strcmp(charges,'2') ...
        & abs(calc_neutral_pepmass-1535.878357)<0.0001 ...
        & num_match_ions >= 8;
    final_result = unique_filtered_result(is_result_selected);

    % sort the results by peptide sequence
    [~,sorted_idx] = sort({final_result.peptide});
    output_res = final_result(sorted_idx);

    % write pep_spec to files
    fileName = fullfile(workspacePath,experimentNames{idx},'pepSpecFile.txt');
    write_peptide_spectra_list_file(output_res, fileName);
end