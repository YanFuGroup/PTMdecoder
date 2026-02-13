function addOrSelectPeptide(obj, sequence)
    % Add a peptide if new, otherwise select the existing peptide
    if isempty(obj.Peptides)
        obj.addPeptide(sequence);
        return;
    end

    existingIdx = find(strcmp({obj.Peptides.peptide_sequence}, sequence), 1, 'first');
    if isempty(existingIdx)
        obj.addPeptide(sequence);
    else
        obj.CurrentPeptideIdx = existingIdx;
        obj.CurrentSpectrumIdx = length(obj.Peptides(existingIdx).spectrum_list);
    end
end
