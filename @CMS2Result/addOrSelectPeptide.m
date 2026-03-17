function addOrSelectPeptide(obj, sequence)
% Add a peptide if new, otherwise select the existing peptide
if isempty(obj.Peptides)
    obj.addPeptide(sequence);
    return;
end

if obj.PeptideIndexMap.Count ~= numel(obj.Peptides)
    obj.rebuildPeptideIndexMap();
end

seq_key = char(sequence);
if isKey(obj.PeptideIndexMap, seq_key)
    existingIdx = double(obj.PeptideIndexMap(seq_key));
    obj.CurrentPeptideIdx = existingIdx;
    obj.CurrentSpectrumIdx = length(obj.Peptides(existingIdx).spectrum_list);
    return;
end

existingIdx = find(strcmp({obj.Peptides.peptide_sequence}, seq_key), 1, 'first');
if isempty(existingIdx)
    obj.addPeptide(sequence);
else
    obj.PeptideIndexMap(seq_key) = int32(existingIdx);
    obj.CurrentPeptideIdx = existingIdx;
    obj.CurrentSpectrumIdx = length(obj.Peptides(existingIdx).spectrum_list);
end
end
