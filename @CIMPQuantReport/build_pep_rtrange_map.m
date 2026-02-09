function pep_rtrange_map = build_pep_rtrange_map(obj)
% Build a map from record key to rt peaks
pep_rtrange_map = containers.Map();
if isempty(obj.blocks)
    return;
end
for idx_block = 1:numel(obj.blocks)
    records = obj.blocks(idx_block).records;
    for idx_rec = 1:numel(records)
        rec = records(idx_rec);
        key = CIMPQuantRecord.build_id(rec.imp_name, rec.charge, rec.raw_name);
        pep_rtrange_map(key) = rec.rt_peaks;
    end
end
end
