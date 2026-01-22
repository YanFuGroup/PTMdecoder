function [isProtN,isProtC] = getWhetherProtNC(obj, strSeq)
% Check whether it is the N-terminus or C-terminus of the protein

% startPos = strfind(obj.m_strOneLine,strSeq);
% if isempty(startPos)
%     % If the peptide sequence cannot be found, report an error and exit
%     error(['Cannot find the peptide: ',strSeq,' in the fasta file!']);
% else
%     % There are two criteria for the N-terminus: one is that the preceding character is the protein termination symbol '$', and the other is that the preceding character is the N-terminal M amino acid of the protein
%     isProtN = obj.m_strOneLine(startPos(1)-1)=='$' || ...
%         (obj.m_strOneLine(startPos(1)-2)=='$'&&obj.m_strOneLine(startPos(1)-1)=='M');
%     isProtC = obj.m_strOneLine(startPos(1)+length(strSeq))=='$';
% end

prot_name_cells = obj.m_pep_prot_mapper.get_proteins(strSeq);
prot_name_pos = obj.get_protein_name_pos(strSeq);

if isempty(prot_name_pos)
    error(['Cannot find the peptide: ', strSeq, ' in the fasta file!']);
else
    % Use the first occurrence to maintain original behavior
    protName = prot_name_pos{1,1};
    pos = prot_name_pos{1,2};
    protSeq = obj.m_mapProt(protName);
    
    % N-terminus criteria: at position 1, or at position 2 with preceding 'M'
    isProtN = (pos == 1) || (pos == 2 && protSeq(1) == 'M');
    % C-terminus criteria: ends at the protein's length
    isProtC = (pos + length(strSeq) - 1) == length(protSeq);
end
end

