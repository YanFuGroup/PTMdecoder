function cfg = buildAlignmentPipelineConfig(obj, aligner, align_strategy, align_options, overrides)
% Build CIMPXICAlignRequantExecutorConfig from current CMSMSPepDeconv fields.
% Input:
%   obj (CMSMSPepDeconv)
%       Processor instance
%   aligner (CXICAligner)
%       Aligner instance
%   align_strategy (CRunAlignStrategy)
%       Alignment strategy
%   align_options (struct, optional)
%       Alignment options
%   overrides (struct, optional)
%       Optional overrides for config fields
% Output:
%   cfg (CIMPXICAlignRequantExecutorConfig)
%       Alignment executor config

if nargin < 4 || isempty(align_options)
    align_options = struct();
end
if nargin < 5 || isempty(overrides)
    overrides = struct();
end

cfg_struct = struct( ...
    'ms12DatasetIO', obj.m_cMs12DatasetIO, ...
    'ms1_tolerance', obj.m_ms1_tolerance, ...
    'minMSMSnum', obj.m_min_MSMS_num, ...
    'alpha', obj.m_alpha, ...
    'resFilterThres', obj.m_resFilterThres, ...
    'aligner', aligner, ...
    'align_strategy', align_strategy, ...
    'align_options', align_options);

override_fields = fieldnames(overrides);
for idx = 1:numel(override_fields)
    cfg_struct.(override_fields{idx}) = overrides.(override_fields{idx});
end

cfg = CIMPXICAlignRequantExecutorConfig(cfg_struct);
end
