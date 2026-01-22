% function cell_prot_name_pos = get_protein_name_pos(obj, pepSeq)
function cell_prot_name_pos = get_protein_name_pos(obj, pepSeq)
% return the protein names containing the specified peptide
% input:
%       pepSeq
%           the sequence of the specified peptide
%       prot_name_cells
%           the name of the protein, cell array of char
% output:
%       cell_prot_name_pos
%           the protein names containing this peptide sequence, and the
%           start position of the peptide on these proteins, in cell
%           format, N*2, N is the number of proteins. [name, pos]

% pos_pepSeq = strfind(obj.m_strOneLine,pepSeq); % find all index
% cell_prot_name_pos = cell(length(pos_pepSeq),2);
% 
% % Parse every protein name according to the position in the single line,
% % find the position of first amino acid is enough.
% for idx = 1:length(pos_pepSeq)
%     [~, parse_idx] = find(pos_pepSeq(idx)>=obj.m_strProtPos,1,'last');
%     cell_prot_name_pos{idx,1} = obj.m_strProtName{parse_idx};
%     cell_prot_name_pos{idx,2} = pos_pepSeq(idx) - obj.m_strProtPos(parse_idx) + 1;
% end


prot_name_cells = obj.m_pep_prot_mapper.get_proteins(pepSeq);
cell_prot_name_pos = {};
for k = 1:length(prot_name_cells)
    prot_name = prot_name_cells{k};
    if isKey(obj.m_mapProt, prot_name)
        seq = obj.m_mapProt(prot_name);
        pos = strfind(seq, pepSeq);
        if ~isempty(pos)
            for i = 1:length(pos)
                cell_prot_name_pos(end+1, :) = {prot_name, pos(i)}; %#ok<AGROW>
            end
        end
    else
        warning(['Protein name ', prot_name, ' not found in the fasta file.']);
    end
end

end