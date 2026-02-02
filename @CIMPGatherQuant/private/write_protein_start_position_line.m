function write_protein_start_position_line(fid, prot_names_pos)
% Write protein start position line.
% Input:
%   fid (1 x 1 double/int)
%       File identifier
%   prot_names_pos (P x 2 cell)
%       Protein name and start position pairs

for idx_np = 1:size(prot_names_pos,1)
    fprintf(fid, '%s,%d;', prot_names_pos{idx_np,1},...
        prot_names_pos{idx_np,2});
end
fprintf(fid,'\n');
end
