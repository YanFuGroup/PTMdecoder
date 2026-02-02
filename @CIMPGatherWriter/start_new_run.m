function start_new_run(output_path)
% start_new_run
% Clear output file and write header for a new run.
% input:
%   output_path (1 x 1 char/string)
%       the path to the output file

fid = fopen(output_path, 'w');
if fid == -1
    error(['Cannot open the the report file ', output_path]);
end
CIMPGatherWriter.write_header(fid);
fclose(fid);
end
