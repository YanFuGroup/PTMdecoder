function raw_names = getRawNamesFromReport(~, quant_report)
% Collect all raw names present in a quantification report.
% Input:
%   quant_report (CIMPQuantReport)
%       Quantification report
% Output:
%   raw_names (1 x N cell)
%       Unique raw names

raw_names = {};
if isempty(quant_report.blocks)
    return;
end

raw_set = containers.Map('KeyType', 'char', 'ValueType', 'logical');
for idx_block = 1:numel(quant_report.blocks)
    records = quant_report.blocks(idx_block).records;
    for idx_rec = 1:numel(records)
        raw_set(records(idx_rec).raw_name) = true;
    end
end
raw_names = raw_set.keys;
end
