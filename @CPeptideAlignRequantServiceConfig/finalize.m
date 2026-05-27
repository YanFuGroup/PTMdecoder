function cfg = finalize(cfg)
% Finalize peptide align-requant config.
% Input:
%   cfg (struct)
%       raw config struct
% Output:
%   cfg (struct)
%       finalized config struct
if ~isfield(cfg, 'align_options') || isempty(cfg.align_options)
    cfg.align_options = struct();
end
if ~isfield(cfg, 'align_strategy_obj') || isempty(cfg.align_strategy_obj)
    CLogger.error(['[CPeptideAlignRequantServiceConfig:InvalidAlignStrategy] ', ...
        'align_strategy_obj must be provided.']);
end
if ~isfield(cfg, 'alignment_report_path')
    cfg.alignment_report_path = '';
end
if ~isfield(cfg, 'requant_output_path')
    cfg.requant_output_path = '';
end
if ~isfield(cfg, 'msms_res_path')
    cfg.msms_res_path = '';
end
if ~isfield(cfg, 'msms_res_multi_file_on') || isempty(cfg.msms_res_multi_file_on)
    cfg.msms_res_multi_file_on = false;
end
if ~isfield(cfg, 'msms_res_paths')
    cfg.msms_res_paths = {};
end
if ~isfield(cfg, 'peptide_quant_res_path')
    cfg.peptide_quant_res_path = '';
end
if ~isfield(cfg, 'align_requant_rt_stats_path')
    cfg.align_requant_rt_stats_path = '';
end
if ~isfield(cfg, 'split_by_dataset_on') || isempty(cfg.split_by_dataset_on)
    cfg.split_by_dataset_on = false;
end
if ~isfield(cfg, 'split_output_dir')
    cfg.split_output_dir = '';
end
if ~isfield(cfg.align_options, 'min_psm') || isempty(cfg.align_options.min_psm)
    cfg.align_options.min_psm = 1;
end
if ~isfield(cfg.align_options, 'num_bins') || isempty(cfg.align_options.num_bins)
    cfg.align_options.num_bins = 5;
end
if ~isfield(cfg.align_options, 'min_per_bin') || isempty(cfg.align_options.min_per_bin)
    cfg.align_options.min_per_bin = 5;
end
if ~isfield(cfg.align_options, 'outlier_k') || isempty(cfg.align_options.outlier_k)
    cfg.align_options.outlier_k = 3;
end
if ~isfield(cfg.align_options, 'outlier_method') || isempty(cfg.align_options.outlier_method)
    cfg.align_options.outlier_method = 'mad';
end
if ~isfield(cfg.align_options, 'rt_sigma') || isempty(cfg.align_options.rt_sigma)
    cfg.align_options.rt_sigma = 0.5;
end
if ~isfield(cfg.align_options, 'max_rt_residual') || isempty(cfg.align_options.max_rt_residual)
    cfg.align_options.max_rt_residual = -1;
end
if ~isfield(cfg.align_options, 'dead_time_min') || isempty(cfg.align_options.dead_time_min)
    cfg.align_options.dead_time_min = 0.5;
end
% Set max_rt_residual to empty if it's negative, using data driven thresholding in that case.
if cfg.align_options.max_rt_residual < 0
    cfg.align_options.max_rt_residual = [];
end

if cfg.msms_res_multi_file_on
    if isempty(cfg.msms_res_paths)
        CLogger.error(['[CPeptideAlignRequantServiceConfig:MissingMsmsResPaths] ', ...
            'msms_res_path_1 ... msms_res_path_N must be provided for peptide align-requant stage.']);
    end
else
    cfg.msms_res_paths = {};
end
if ~cfg.msms_res_multi_file_on && isempty(cfg.msms_res_path)
    CLogger.error(['[CPeptideAlignRequantServiceConfig:MissingMsmsResPath] ', ...
        'msms_res_path must be provided for peptide align-requant stage.']);
end
if isempty(cfg.peptide_quant_res_path)
    CLogger.error(['[CPeptideAlignRequantServiceConfig:MissingPeptideQuantResPath] ', ...
        'peptide_quant_res_path must be provided for peptide align-requant stage.']);
end
if isempty(cfg.alignment_report_path)
    CLogger.error(['[CPeptideAlignRequantServiceConfig:MissingAlignReportPath] ', ...
        'align_report_path must be provided for peptide align-requant stage.']);
end
if isempty(cfg.requant_output_path)
    CLogger.error(['[CPeptideAlignRequantServiceConfig:MissingAlignRequantOutputPath] ', ...
        'align_requant_output_path must be provided for peptide align-requant stage.']);
end
if cfg.split_by_dataset_on && isempty(cfg.split_output_dir)
    combined_output_dir = fileparts(char(string(cfg.requant_output_path)));
    cfg.split_output_dir = fullfile(combined_output_dir, 'split_by_dataset');
end
end
