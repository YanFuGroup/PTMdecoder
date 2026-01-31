function [isProtN,isProtC] = getWhetherProtNC(obj, strSeq)
% Check whether it is the N-terminus or C-terminus of the protein
% Input:
%   obj (CPepProtService)
%       Protein/peptide service instance
%   strSeq (1 x 1 char/string)
%       Peptide sequence
% Output:
%   isProtN (1 x 1 logical)
%       Whether peptide is at protein N-terminus
%   isProtC (1 x 1 logical)
%       Whether peptide is at protein C-terminus

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
