clear
% syn_ind shared preprocess pepSpec selection matches preprocess_syn_mix.m:
% +2 precursor charge and >= 8 Mascot matched fragment ions.
%% User configuration
% Update these paths before running this script.
% The input and output folders do not need to share the same parent folder.
res_path = '<path_to_syn_ind_mascot_dat>';
work_dir_res = '<path_to_syn_ind_preprocess_result>';
experimentNames = {'ind1','ind2','ind3','ind4','ind5','ind6','ind7','ind8','ind9','ind10', ...
    'ind11','ind12','ind13','ind14','ind15','ind16'};

if strcmp(res_path,'<path_to_syn_ind_mascot_dat>') || ~isfolder(res_path)
    error('Please set res_path to an existing Mascot DAT folder before running this script: %s', res_path);
end

if strcmp(work_dir_res,'<path_to_syn_ind_preprocess_result>')
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

    % select the target peptides GKGGKGLGKGGAKR
    peptides = {unique_filtered_result.peptide}';
    charge_values = str2double({unique_filtered_result.Charge}');
    calc_neutral_pepmass = [unique_filtered_result.Calc_neutral_pepmass]';
    num_match_ions = str2double({unique_filtered_result.num_match_ions}');
    is_result_selected = strcmp(peptides,'GKGGKGLGKGGAKR') ...
        & charge_values == 2 ...
        & abs(calc_neutral_pepmass-1535.878357)<0.0001 ...
        & num_match_ions >= 8;
    final_result = unique_filtered_result(is_result_selected);

    % sort the results by peptide sequence
    [~,sorted_idx] = sort({final_result.peptide});
    output_res = final_result(sorted_idx);

    % write pep_spec to files
    fileName = fullfile(outputDir,'pepSpecFile.txt');
    write_peptide_spectra_list_file(output_res, fileName);
end
