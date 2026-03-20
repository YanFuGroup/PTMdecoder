function write_protein_start_position_line(fid, prot_names_pos)
% write_protein_start_position_line
% Write protein name and start position line.
% input:
%   fid (1 x 1 double)
%       file identifier for the output file
%   prot_names_pos (N x 2 cell)
%       cell array containing protein names and their start positions

if isempty(prot_names_pos)
    fprintf(fid, '\n');
    return;
end

% for idx_np = 1:size(prot_names_pos,1)
%     fprintf(fid, '%s,%d;', prot_names_pos{idx_np,1}, ...
%         prot_names_pos{idx_np,2});
% end
% fprintf(fid,'\n');

prot_names_pos_T = prot_names_pos';
str = sprintf('%s,%d;', prot_names_pos_T{:});
fprintf(fid, '%s\n', str);
end
