function write(result, path)
% write - Minimal writer for MS2 processing result struct

fid = fopen(path, 'w');
if fid <= 0
    error(['Cannot open output file: ', path]);
end
cleanup = onCleanup(@() fclose(fid));

if isfield(result, 'header') && ~isempty(result.header)
    fprintf(fid, '%s\n', result.header);
end

if isfield(result, 'lines') && ~isempty(result.lines)
    for i = 1:numel(result.lines)
        fprintf(fid, '%s\n', result.lines{i});
    end
end
end
