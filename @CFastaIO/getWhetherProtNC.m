function [isProtN,isProtC] = getWhetherProtNC(obj,strSeq)
% Check whether it is the N-terminus or C-terminus of the protein
startPos = strfind(obj.m_strOneLine,strSeq);
if isempty(startPos)
    % If the peptide sequence cannot be found, report an error and exit
    error(['Cannot find the peptide: ',strSeq,' in the fasta file!']);
else
    % There are two criteria for the N-terminus: one is that the preceding character is the protein termination symbol '$', and the other is that the preceding character is the N-terminal M amino acid of the protein
    isProtN = obj.m_strOneLine(startPos(1)-1)=='$' || ...
        (obj.m_strOneLine(startPos(1)-2)=='$'&&obj.m_strOneLine(startPos(1)-1)=='M');
    isProtC = obj.m_strOneLine(startPos(1)+length(strSeq))=='$';
end
end

