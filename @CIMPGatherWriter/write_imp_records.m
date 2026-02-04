function write_imp_records(fid, imp_records)
% write_imp_records
% Write IMP records to an open file handle.
% input:
%   fid (1 x 1 double)
%       file identifier for the output file
%   imp_records (struct array or CIMPQuantResult array)
%       array of structures containing IMP record data

if isempty(imp_records)
    return;
end

for idx_rec = 1:numel(imp_records)
    rec = imp_records(idx_rec);
    imp_name = rec.imp_name;
    raw_name = rec.raw_name;
    fprintf(fid, '*\t%s\t+%d\t%s\t%.4f\t%f\t%f\t%f\n', ...
        imp_name, rec.charge, raw_name, rec.mass_center, ...
        rec.low_mz_bound, rec.high_mz_bound, rec.area);
    for idx_peak = 1:numel(rec.rt_peaks)
        peak = rec.rt_peaks(idx_peak);
        fprintf(fid, '@\t%f\t%f\t%f\t%d\n', ...
            peak.rt_start, peak.rt_end, peak.ratio, peak.check_label);
    end
end
end
