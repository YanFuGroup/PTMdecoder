function write(result, path)
% write - Write MS2 processing result to file
% Input:
%   result (CMS2Result)
%       Write hierarchical result as P/S/peptidoform lines with named fields.
%   path (1 x N char/string)
%       Output file path

fieldKeys = CMS2ResultIO.CMS2ResultFieldKeys();

fid = fopen(path, 'w');
if fid <= 0
    error('CMS2ResultIO:OpenFileFailed', 'Cannot open output file: %s', path);
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
        jaccardVal = NaN;
        if isfield(spectrum, 'jaccard_stability') && ~isempty(spectrum.jaccard_stability)
            jaccardVal = spectrum.jaccard_stability;
        end

        fprintf(fid, 'S\t%s\t%s\t%s=%.6f\n', ...
            spectrum.dataset_name, ...
            spectrum.spectrum_name, ...
            fieldKeys.spectrum.jaccard, ...
            jaccardVal);

        for k = 1:spectrum.peptidoform_num
            supportVal = NaN;
            madVal = NaN;
            if isfield(spectrum, 'peptidoform_list_support_freq') && ...
                    numel(spectrum.peptidoform_list_support_freq) >= k && ...
                    ~isempty(spectrum.peptidoform_list_support_freq(k))
                supportVal = spectrum.peptidoform_list_support_freq(k);
            end
            if isfield(spectrum, 'peptidoform_list_abundance_mad') && ...
                    numel(spectrum.peptidoform_list_abundance_mad) >= k && ...
                    ~isempty(spectrum.peptidoform_list_abundance_mad(k))
                madVal = spectrum.peptidoform_list_abundance_mad(k);
            end

            fprintf(fid, '%s\t%.6f\t%s=%.6f\t%s=%.6f\n', ...
                spectrum.peptidoform_list_str{k}, ...
                spectrum.peptidoform_list_abun(k), ...
                fieldKeys.peptidoform.support, ...
                supportVal, ...
                fieldKeys.peptidoform.mad, ...
                madVal);
        end
    end
end
end
