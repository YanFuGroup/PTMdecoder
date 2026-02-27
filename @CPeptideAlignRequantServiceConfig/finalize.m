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
if ~isfield(cfg, 'alignment_report_path')
    cfg.alignment_report_path = '';
end
if ~isfield(cfg, 'requant_output_path')
    cfg.requant_output_path = '';
end
% Set max_rt_residual to empty if it's negative, using data driven thresholding in that case.
if cfg.align_options.max_rt_residual < 0
    cfg.align_options.max_rt_residual = [];
end

cfg.alignment_report_path = CPathResolver.resolveFilePath( ...
    cfg.output_dir_path, 'report_alignment.txt', cfg.alignment_report_path);
cfg.requant_output_path = CPathResolver.resolveFilePath( ...
    cfg.output_dir_path, 'report_peptide_all_requant_aligned.txt', cfg.requant_output_path);
end