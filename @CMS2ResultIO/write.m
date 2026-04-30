function write(result, path, include_vif)
% write - Write MS2 processing result to file
% Input:
%   result (CMS2Result)
%       Write hierarchical result as P/S/peptidoform lines with named fields.
%   path (1 x N char/string)
%       Output file path
%   include_vif (1 x 1 logical, optional)
%       True to include VIF named fields in the output report.

if nargin < 3 || isempty(include_vif)
    include_vif = false;
end
include_vif = logical(include_vif);

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
        vifAllVal = NaN;
        vifReportedVal = NaN;
        if isfield(spectrum, 'jaccard_stability') && ~isempty(spectrum.jaccard_stability)
            jaccardVal = spectrum.jaccard_stability;
        end
        if isfield(spectrum, 'vif_all_imp_max') && ~isempty(spectrum.vif_all_imp_max)
            vifAllVal = spectrum.vif_all_imp_max;
        end
        if isfield(spectrum, 'vif_reported_imp_max') && ~isempty(spectrum.vif_reported_imp_max)
            vifReportedVal = spectrum.vif_reported_imp_max;
        end

        if include_vif
            fprintf(fid, 'S\t%s\t%s\t%s=%.6f\t%s=%.6f\t%s=%.6f\n', ...
                spectrum.dataset_name, ...
                spectrum.spectrum_name, ...
                fieldKeys.spectrum.jaccard, ...
                jaccardVal, ...
                fieldKeys.spectrum.vifAll, ...
                vifAllVal, ...
                fieldKeys.spectrum.vifReported, ...
                vifReportedVal);
        else
            fprintf(fid, 'S\t%s\t%s\t%s=%.6f\n', ...
                spectrum.dataset_name, ...
                spectrum.spectrum_name, ...
                fieldKeys.spectrum.jaccard, ...
                jaccardVal);
        end

        for k = 1:spectrum.peptidoform_num
            supportVal = NaN;
            vifVal = NaN;
            madVal = NaN;
            if isfield(spectrum, 'peptidoform_list_support_freq') && ...
                    numel(spectrum.peptidoform_list_support_freq) >= k && ...
                    ~isempty(spectrum.peptidoform_list_support_freq(k))
                supportVal = spectrum.peptidoform_list_support_freq(k);
            end
            if isfield(spectrum, 'peptidoform_list_vif') && ...
                    numel(spectrum.peptidoform_list_vif) >= k && ...
                    ~isempty(spectrum.peptidoform_list_vif(k))
                vifVal = spectrum.peptidoform_list_vif(k);
            end
            if isfield(spectrum, 'peptidoform_list_abundance_mad') && ...
                    numel(spectrum.peptidoform_list_abundance_mad) >= k && ...
                    ~isempty(spectrum.peptidoform_list_abundance_mad(k))
                madVal = spectrum.peptidoform_list_abundance_mad(k);
            end

            if include_vif
                fprintf(fid, '%s\t%.6f\t%s=%.6f\t%s=%.6f\t%s=%.6f\n', ...
                    spectrum.peptidoform_list_str{k}, ...
                    spectrum.peptidoform_list_abun(k), ...
                    fieldKeys.peptidoform.support, ...
                    supportVal, ...
                    fieldKeys.peptidoform.vif, ...
                    vifVal, ...
                    fieldKeys.peptidoform.mad, ...
                    madVal);
            else
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
end
