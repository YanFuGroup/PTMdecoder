function write_mascot_result_table(result, filename)
% Write Mascot result table to file
% Input:
%   result (1 x N struct)
%       Mascot result entries
%   filename (1 x 1 char/string)
%       output file path
CFdrFilteredResultIO.write(result, filename);
end