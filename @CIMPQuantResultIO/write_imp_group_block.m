function write_imp_group_block(fid, protein_name_pos, imp_records)
% write_imp_group_block
% Append IMP group block to an existing output file.
% input:
%   fid (1 x 1 double)
%       file identifier for the output file
%   protein_name_pos (N x 2 cell)
%       cell array containing protein names and their start positions
%   imp_records (struct array or CIMPQuantRecord array)
%       array of structures containing IMP record data

if isempty(imp_records)
    return;
end
CIMPQuantResultIO.write_protein_start_position_line(fid, protein_name_pos);
CIMPQuantResultIO.write_imp_records(fid, imp_records);
end
