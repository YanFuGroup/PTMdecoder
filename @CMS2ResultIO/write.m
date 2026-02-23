function write(result, path)
% write - Write MS2 processing result to file
% Input:
%   result (CMS2Result)
%       Write hierarchical result as P/S/peptidoform lines.
%   path (1 x N char/string)
%       Output file path

fid = fopen(path, 'w');
if fid <= 0
    error(['Cannot open output file: ', path]);
end
cleanup = onCleanup(@() fclose(fid));

if ~isa(result, 'CMS2Result')
    error('CMS2ResultIO:InvalidWriteInput', ...
        'write expects a CMS2Result object.');
end

for i = 1:length(result.Peptides)
    if i > 1
        fprintf(fid, '\n\n');
    end
    peptide = result.Peptides(i);
    fprintf(fid, 'P\t%s\n', peptide.peptide_sequence);
    for j = 1:length(peptide.spectrum_list)
        spectrum = peptide.spectrum_list(j);
        fprintf(fid, 'S\t%s\t%s\n', spectrum.dataset_name, spectrum.spectrum_name);
        for k = 1:spectrum.peptidoform_num
            fprintf(fid, '%s\t%.6f\n', spectrum.peptidoform_list_str{k}, spectrum.peptidoform_list_abun(k));
        end
    end
end
end
