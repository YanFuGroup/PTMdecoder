function log_quant_group_summary(stage_name)
% Write an INFO summary for initial quantification group outcomes.

stats = CIMPQuantStats.quant_group_stats('get', []);
field_names = fieldnames(stats);
parts = cell(size(field_names));
for idx = 1:numel(field_names)
    parts{idx} = sprintf('%s=%d', field_names{idx}, stats.(field_names{idx}));
end
CLogger.info('%s quant group summary: %s.', stage_name, strjoin(parts, ', '));
end
