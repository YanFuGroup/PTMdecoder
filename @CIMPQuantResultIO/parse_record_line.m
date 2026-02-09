function rec = parse_record_line(strline)
% Parse a record line starting with '*'
% Line format: *\tIMP_name\t+charge\tdataset_name\tmass_center\tlow_mz_bound\thigh_mz_bound\tarea
% Exactly 8 tab-separated fields, starting with '*'
segment = regexp(strline, '\t', 'split');
imp_name = segment{2};
charge = str2double(strrep(segment{3}, '+', ''));
raw_name = segment{4};
mass_center = str2double(segment{5});
low_mz_bound = str2double(segment{6});
high_mz_bound = str2double(segment{7});
area = str2double(segment{8});
rt_peaks = repmat(struct('rt_start',0,'rt_end',0,'ratio',0,'check_label',0), 0, 1);
rec = CIMPQuantRecord(imp_name, charge, raw_name, mass_center, low_mz_bound, high_mz_bound, area, rt_peaks);
end
