function write(report, path)
% Write report content to file (format compatible with legacy output)
fid = fopen(path, 'w');
if fid < 0
    CLogger.error(['Cannot open the report file "', path, '".']);
end
cleanup_output = onCleanup(@() fclose(fid));
CIMPQuantResultIO.write_header(fid);

if isempty(report.blocks)
    return;
end
for idx_block = 1:numel(report.blocks)
    block = report.blocks(idx_block);
    CIMPQuantResultIO.write_imp_group_block(fid, block.protein_name_pos, block.records);
end
clear cleanup_output;
end
