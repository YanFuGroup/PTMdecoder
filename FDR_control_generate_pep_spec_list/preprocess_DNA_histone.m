clear
%% User configuration
% Update these paths before running this script.
% The input and output folders do not need to share the same parent folder.
res_path = '<path_to_DNA_histone_mascot_dat>';
work_dir_res = '<path_to_DNA_histone_results>';
experimentNames = {'wt1','wt2','wt3','ir1','ir2','ir3'};

if strcmp(res_path,'<path_to_DNA_histone_mascot_dat>') || ~isfolder(res_path)
    error('Please set res_path to an existing Mascot DAT folder before running this script: %s', res_path);
end

if strcmp(work_dir_res,'<path_to_DNA_histone_results>')
    error('Please set work_dir_res before running this script.');
end

if ~isfolder(work_dir_res)
    [status,message] = mkdir(work_dir_res);
    if ~status
        error('Failed to create output folder "%s": %s', work_dir_res, message);
    end
end

for idx = 1:length(experimentNames)
    outputDir = fullfile(work_dir_res,experimentNames{idx});
    if ~isfolder(outputDir)
        [status,message] = mkdir(outputDir);
        if ~status
            error('Failed to create experiment output folder "%s": %s', outputDir, message);
        end
    end

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
    grouped_filename = fullfile(outputDir,'group_result_mascot.txt');
    write_mascot_result_table(group_filtered_result, grouped_filename);

    % Write all filtered result
    FDR_filtered_result = result(Iid.SF);
    filtered_fileName = fullfile(outputDir,'filtered_result_mascot.txt');
    write_mascot_result_table(FDR_filtered_result, filtered_fileName);

    % select the unique identified spectra
    mgf_scan_strings = strcat({FDR_filtered_result.DatasetName}',' ',{FDR_filtered_result.Scan}');
    [~,~,group_idx] = unique(mgf_scan_strings,'stable');
    group_counts = accumarray(group_idx,1);
    is_result_selected = (group_counts(group_idx) == 1);
    unique_filtered_result = FDR_filtered_result(is_result_selected);

    % Select the modified peptides without Oxidation, cannot be all fixed modification
    modifications = {unique_filtered_result.modification}';
    is_not_dash = ~strcmp(modifications,'-');
    is_not_oxidation = ~contains(modifications,'Oxidation');
    candidate_idx = is_not_dash & is_not_oxidation;
    is_all_label_only = false(size(modifications));
    is_all_label_only(candidate_idx) = cellfun(@(mod_text) all(ismember(strsplit(mod_text, ','), ...
        {'Label:13C(6) (K)','Label:13C(6)15N(4) (R)'})), modifications(candidate_idx));
    is_result_selected = candidate_idx & ~is_all_label_only;
    final_result = unique_filtered_result(is_result_selected);

    % sort the results by peptide sequence
    [~,sorted_idx] = sort({final_result.peptide});
    output_res = final_result(sorted_idx);

    % write pep_spec to files
    fileName = fullfile(outputDir,'pepSpecFile.txt');
    write_peptide_spectra_list_file(output_res, fileName);

end
