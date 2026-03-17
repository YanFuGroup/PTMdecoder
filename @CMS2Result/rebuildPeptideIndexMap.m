function rebuildPeptideIndexMap(obj)
% Rebuild peptide-sequence to index cache from current Peptides array
obj.PeptideIndexMap = containers.Map('KeyType', 'char', 'ValueType', 'int32');
for idxPep = 1:numel(obj.Peptides)
    seq_key = char(obj.Peptides(idxPep).peptide_sequence);
    obj.PeptideIndexMap(seq_key) = int32(idxPep);
end
end
