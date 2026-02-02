function write_imp_group_block(output_path, protein_name_pos, imp_records)
% write_imp_group_block
% Append IMP group block to an existing output file.
% input:
%   output_path (1 x 1 char/string)
%       the path to the output file
%   protein_name_pos (N x 2 cell)
%       cell array containing protein names and their start positions
%   imp_records (struct array)
%       array of structures containing IMP record data

if isempty(imp_records)
    return;
end
fid = fopen(output_path, 'a');
if fid == -1
    error(['Cannot open the the report file ', output_path]);
end
CIMPGatherWriter.write_protein_start_position_line(fid, protein_name_pos);
CIMPGatherWriter.write_imp_records(fid, imp_records);
fclose(fid);
end
