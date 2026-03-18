function cfg = fromParamMap(task_param_map)
% Build site-level pipeline config from parameter map.
% Input:
%   task_param_map (containers.Map)
%       parameter key-value map for site-level summary
% Output:
%   cfg (CSiteLevelPipelineConfig)
if ~isa(task_param_map, 'containers.Map')
    error('CSiteLevelPipelineConfig:InvalidParamMap', ...
        'Expected task_param_map to be a containers.Map.');
end

output_dir_path = CParamMapUtils.getOptional(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_OUTPUT_DIR_PATH, '', 'CSiteLevelPipelineConfig');
if ~isempty(output_dir_path)
    pep_default = fullfile(output_dir_path, 'report_peptide_all.txt');
    intere_default = fullfile(output_dir_path, 'report_site.txt');
    unintere_default = fullfile(output_dir_path, 'report_peptide_uninterested.txt');
    if ~exist(output_dir_path, 'dir')
        mkdir(output_dir_path);
    end
else
    pep_default = '';
    intere_default = '';
    unintere_default = '';
end

cfg_struct = struct();
cfg_struct.input_path = CParamMapUtils.getOptional(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_PEP_LEVEL_FILE_PATH, pep_default, 'CSiteLevelPipelineConfig');
cfg_struct.output_intere_path = CParamMapUtils.getOptional(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_OUTPUT_INTERE_PATH, intere_default, 'CSiteLevelPipelineConfig');
cfg_struct.output_unintere_path = CParamMapUtils.getOptional(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_OUTPUT_UNINTERE_PATH, unintere_default, 'CSiteLevelPipelineConfig');
protein_abbr_input_mode = CParamMapUtils.getOptional( ...
    task_param_map, ...
    CPTMdecoderWorkflowParamKeys.PARAM_SITE_LEVEL_PROTEIN_ABBR_INPUT_MODE, ...
    'manual', ...
    'CSiteLevelPipelineConfig');
if isstring(protein_abbr_input_mode)
    protein_abbr_input_mode = char(protein_abbr_input_mode);
end
protein_abbr_input_mode = lower(strtrim(protein_abbr_input_mode));
cfg_struct.protein_abbr_input_mode = protein_abbr_input_mode;

switch protein_abbr_input_mode
    case 'manual'
        protein_name_abbr_num = CParamMapUtils.getRequiredNumber( ...
            task_param_map, ...
            CPTMdecoderWorkflowParamKeys.PARAM_PROTEIN_NAME_ABBR_NUM, ...
            'number of protein abbreviation mappings');
        protein_name_abbr = containers.Map;
        for idx = 1:protein_name_abbr_num
            key_name = [CPTMdecoderWorkflowParamKeys.PARAM_PREFIX_PROTEIN_NAME_ABBR, num2str(idx)];
            pair_str = CParamMapUtils.getRequired(task_param_map, key_name, 'protein abbreviation pair');
            split_str = strsplit(pair_str, '>');
            protein_name_abbr(strtrim(split_str{1})) = strtrim(split_str{2});
        end

        cfg_struct.protein_name_abbr = protein_name_abbr;

    case 'file'
        protein_abbr_file_path = CParamMapUtils.getRequired( ...
            task_param_map, ...
            CPTMdecoderWorkflowParamKeys.PARAM_SITE_LEVEL_PROTEIN_ABBR_FILE_PATH, ...
            'site-level protein abbreviation map file path', ...
            'CSiteLevelPipelineConfig');
        protein_abbr_file_col_protein_name = CParamMapUtils.getRequired( ...
            task_param_map, ...
            CPTMdecoderWorkflowParamKeys.PARAM_SITE_LEVEL_PROTEIN_ABBR_FILE_COL_PROTEIN_NAME, ...
            'site-level protein abbreviation file protein-name column', ...
            'CSiteLevelPipelineConfig');
        protein_abbr_file_col_abbr_name = CParamMapUtils.getRequired( ...
            task_param_map, ...
            CPTMdecoderWorkflowParamKeys.PARAM_SITE_LEVEL_PROTEIN_ABBR_FILE_COL_ABBR_NAME, ...
            'site-level protein abbreviation file abbreviation column', ...
            'CSiteLevelPipelineConfig');

        cfg_struct.protein_name_abbr = CSiteLevelPipelineConfig.loadProteinAbbrMapFromTsv( ...
            protein_abbr_file_path, ...
            protein_abbr_file_col_protein_name, ...
            protein_abbr_file_col_abbr_name);
        cfg_struct.protein_abbr_file_path = protein_abbr_file_path;
        cfg_struct.protein_abbr_file_col_protein_name = protein_abbr_file_col_protein_name;
        cfg_struct.protein_abbr_file_col_abbr_name = protein_abbr_file_col_abbr_name;

    otherwise
        error('CSiteLevelPipelineConfig:InvalidProteinAbbrInputMode', ...
            'site-level protein abbreviation input mode must be either ''manual'' or ''file'' (got ''%s'').', ...
            protein_abbr_input_mode);
end

mod_name_abbr_num = CParamMapUtils.getRequiredNumber(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_MOD_NAME_ABBR_NUM, 'number of modification abbreviation mappings');
mod_name_abbr = containers.Map;
for idx = 1:mod_name_abbr_num
    key_name = [CPTMdecoderWorkflowParamKeys.PARAM_PREFIX_MOD_NAME_ABBR, num2str(idx)];
    pair_str = CParamMapUtils.getRequired(task_param_map, key_name, 'modification abbreviation pair');
    split_str = strsplit(pair_str, '>');
    mod_name_abbr(strtrim(split_str{1})) = strtrim(split_str{2});
end
cfg_struct.mod_name_abbr = mod_name_abbr;

ignore_str = CParamMapUtils.getRequired(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_IGNORE_STRINGS_SITE_LEVEL, 'ignore strings for site-level summary');
cfg_struct.ignore_strings = CParamMapUtils.parseQuotedList(ignore_str, ';');
cfg_struct.column_idxs = struct('icol_seq', 2, 'icol_auc', 8);

if isempty(cfg_struct.input_path)
    error('CSiteLevelPipelineConfig:InvalidInputPath', 'input_path must be provided.');
end
if isempty(cfg_struct.output_intere_path)
    error('CSiteLevelPipelineConfig:InvalidInterestedOutputPath', 'output_intere_path must be provided.');
end
if isempty(cfg_struct.output_unintere_path)
    error('CSiteLevelPipelineConfig:InvalidUninterestedOutputPath', 'output_unintere_path must be provided.');
end
if ~isfield(cfg_struct.column_idxs, 'icol_seq') || ~isfield(cfg_struct.column_idxs, 'icol_auc')
    error('CSiteLevelPipelineConfig:InvalidColumnIdxs', 'column_idxs must have fields icol_seq and icol_auc.');
end

cfg = CSiteLevelPipelineConfig(cfg_struct);
end
