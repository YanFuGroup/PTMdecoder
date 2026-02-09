function pep_prot_map = build_pep_prot_map(obj)
% Build a map from record key to the first protein name in the block
pep_prot_map = containers.Map();
if isempty(obj.blocks)
    return;
end
for idx_block = 1:numel(obj.blocks)
    prot_name = '';
    if ~isempty(obj.blocks(idx_block).protein_name_pos)
        prot_name = obj.blocks(idx_block).protein_name_pos{1,1};
    end
    records = obj.blocks(idx_block).records;
    for idx_rec = 1:numel(records)
        rec = records(idx_rec);
        key = CIMPQuantRecord.build_id(rec.imp_name, rec.charge, rec.raw_name);
        pep_prot_map(key) = prot_name;
    end
end
end
