function cell_prot_name_pos = get_protein_name_pos(obj, pepSeq)
% return the protein names containing the specified peptide
% input:
%       pepSeq
%           the sequence of the specified peptide
% output:
%       cell_prot_name_pos
%           the protein names containing this peptide sequence, and the
%           start position of the peptide on these proteins, in cell
%           format, N*2, N is the number of proteins. [name, pos]

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
