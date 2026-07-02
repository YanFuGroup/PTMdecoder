clear
%% User configuration
% Update these paths before running this script.
% The input and output folders do not need to share the same parent folder.
res_path = '<path_to_MCF7_JIB04_mascot_dat>';
work_dir_res = '<path_to_MCF7_JIB04_result>';
work_dir_res_top1 = '<path_to_MCF7_JIB04_result_top1>';
experimentNames = {'DMSO_1','DMSO_2','JIB04_1','JIB04_2'};

if strcmp(res_path,'<path_to_MCF7_JIB04_mascot_dat>') || ~isfolder(res_path)
    error('Please set res_path to an existing Mascot DAT folder before running this script: %s', res_path);
end

if strcmp(work_dir_res,'<path_to_MCF7_JIB04_result>')
    error('Please set work_dir_res before running this script.');
end

if strcmp(work_dir_res_top1,'<path_to_MCF7_JIB04_result_top1>')
    error('Please set work_dir_res_top1 before running this script.');
end

if ~isfolder(work_dir_res)
    [status,message] = mkdir(work_dir_res);
    if ~status
        error('Failed to create output folder "%s": %s', work_dir_res, message);
    end
end

if ~isfolder(work_dir_res_top1)
    [status,message] = mkdir(work_dir_res_top1);
    if ~status
        error('Failed to create output folder "%s": %s', work_dir_res_top1, message);
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

    outputDirTop1 = fullfile(work_dir_res_top1,experimentNames{idx});
    if ~isfolder(outputDirTop1)
        [status,message] = mkdir(outputDirTop1);
        if ~status
            error('Failed to create experiment output folder "%s": %s', outputDirTop1, message);
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
    grouped_filename = fullfile(outputDirTop1,'group_result_mascot.txt');
    write_mascot_result_table(group_filtered_result, grouped_filename);

    % Write all filtered result
    FDR_filtered_result = result(Iid.SF);
    filtered_fileName = fullfile(outputDir,'filtered_result_mascot.txt');
    write_mascot_result_table(FDR_filtered_result, filtered_fileName);
    filtered_fileName = fullfile(outputDirTop1,'filtered_result_mascot.txt');
    write_mascot_result_table(FDR_filtered_result, filtered_fileName);

    % select the unique identified spectra
    mgf_scan_strings = strcat({FDR_filtered_result.DatasetName}',' ',{FDR_filtered_result.Scan}');
    [~,~,group_idx] = unique(mgf_scan_strings,'stable');
    group_counts = accumarray(group_idx,1);
    is_result_selected = (group_counts(group_idx) == 1);
    unique_filtered_result = FDR_filtered_result(is_result_selected);

    % select the modified peptides without Oxidation
    modifications = {unique_filtered_result.modification}';
    is_result_selected = ~strcmp(modifications,'-') & ~contains(modifications,'Oxidation');
    final_result = unique_filtered_result(is_result_selected);

    % sort the results by peptide sequence
    [~,sorted_idx] = sort({final_result.peptide});
    output_res = final_result(sorted_idx);

    % write pep_spec to files
    fileName = fullfile(outputDir,'pepSpecFile.txt');
    write_peptide_spectra_list_file(output_res, fileName);
    % write the top scoring results to file in report_msms.txt format
    fileName = fullfile(outputDirTop1,'report_msms_top1.txt');
    write_report_msms_top1(output_res, fileName);

end
