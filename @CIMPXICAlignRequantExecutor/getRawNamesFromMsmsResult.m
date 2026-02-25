function raw_names = getRawNamesFromMsmsResult(~, msms_result)
% Collect raw names from MSMS results.
% Input:
%   msms_result (CMS2Result)
%       MSMS results from report_msms.txt
% Output:
%   raw_names (1 x N cell)
%       Unique raw names

raw_set = containers.Map('KeyType', 'char', 'ValueType', 'logical');
for idx_pep = 1:length(msms_result.Peptides)
    spectrum_list = msms_result.Peptides(idx_pep).spectrum_list;
    for idx_spec = 1:length(spectrum_list)
        raw_set(spectrum_list(idx_spec).dataset_name) = true;
    end
end
raw_names = raw_set.keys;
end
