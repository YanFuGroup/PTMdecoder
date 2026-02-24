function cfg = buildIMPProcessingPipelineConfig(obj, overrides)
% Build CIMPProcessingPipelineConfig from current CMSMSPepDeconv fields.
% Input:
%   obj (CMSMSPepDeconv)
%       Processor instance
%   overrides (struct, optional)
%       Optional overrides for fields:
%       - ms12DatasetIO
%       - ms1_tolerance
%       - minMSMSnum
%       - alpha
%       - resFilterThres
% Output:
%   cfg (CIMPProcessingPipelineConfig)
%       Processing pipeline config

if nargin < 2 || isempty(overrides)
    overrides = struct();
end

cfg_struct = struct( ...
    'ms12DatasetIO', obj.m_cMs12DatasetIO, ...
    'ms1_tolerance', obj.m_ms1_tolerance, ...
    'minMSMSnum', obj.m_min_MSMS_num, ...
    'alpha', obj.m_alpha, ...
    'resFilterThres', obj.m_resFilterThres);

override_fields = fieldnames(overrides);
for idx = 1:numel(override_fields)
    cfg_struct.(override_fields{idx}) = overrides.(override_fields{idx});
end

cfg = CIMPProcessingPipelineConfig(cfg_struct);
end
