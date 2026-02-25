function raw_names = getRawNamesFromRawIdentManagers(~, rawIdentManagers)
% Collect raw names from prebuilt raw identification managers.
% Input:
%   rawIdentManagers (1 x N cell)
%       Per-peptide CIMPRawIdentManager list
% Output:
%   raw_names (1 x N cell)
%       Unique raw names

raw_set = containers.Map('KeyType', 'char', 'ValueType', 'logical');
for idx_pep = 1:length(rawIdentManagers)
    rawIdentManager = rawIdentManagers{idx_pep};
    if isempty(rawIdentManager)
        continue;
    end
    if ~isa(rawIdentManager, 'CIMPRawIdentManager')
        error('getRawNamesFromRawIdentManagers:InvalidInput', ...
            'rawIdentManagers{%d} must be CIMPRawIdentManager.', idx_pep);
    end
    cur_names = rawIdentManager.getEntries();
    for idx_raw = 1:length(cur_names)
        raw_set(cur_names{idx_raw}) = true;
    end
end
raw_names = raw_set.keys;
end
