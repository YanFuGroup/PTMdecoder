function write(filteredResults, filename)
% Write a result table to file.
% Input:
%   filteredResults (CFdrFilteredResult or 1 x N struct)
%       Entries to write. If a container is provided, uses its entries.
%   filename (1 x 1 char/string)
%       Output file path.
%
% The output is a tab-delimited table with a header line and columns in the
% same order as the FDR filtered result table.

if isa(filteredResults, 'CFdrFilteredResult')
    entries = filteredResults.entries;
else
    entries = filteredResults;
end

fout = fopen(filename, 'w');
if fout == -1
    error('Failed to open the file "%s"!', filename);
end
cleanup = onCleanup(@() fclose(fout));

fprintf(fout, ['Site\tDatasetName\tScan\tSpectrum\tCharge\t' ...
    'Calc_neutral_pepmass\tprecursor_neutral_mass\tmassdiff' ...
    '\tnum_match_ions\tpeptide\tprotein\tmodification\t' ...
    'modificationlocation\tScore\n']);

for i_entry = 1:length(entries)
    e = entries(i_entry);
    fields = { ...
        e.Site, e.DatasetName, e.Scan, e.Spectrum, e.Charge, ...
        e.Calc_neutral_pepmass, e.precursor_neutral_mass, e.massdiff, ...
        e.num_match_ions, e.peptide, e.protein, e.modification, ...
        e.modificationlocation, e.Score ...
    };
    fields = cellfun(@CFdrFilteredResultIO.toString, fields, 'UniformOutput', false);
    fprintf(fout, '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n', fields{:});
end
end
