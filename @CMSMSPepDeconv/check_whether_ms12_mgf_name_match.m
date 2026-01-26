function check_whether_ms12_mgf_name_match(obj)
% Check whether the name of ms1/ms2 is matched with mgf

mgf_files = dir(fullfile(obj.m_specPath, '*.mgf'));
ms1_files = dir(fullfile(obj.m_specPath, '*.ms1'));

% Get the stem name of the files
mgf_names = {mgf_files.name};
ms1_names = {ms1_files.name};

% Remove the extension
mgf_names = cellfun(@(x) x(1:end-4), mgf_names, 'UniformOutput', false);
ms1_names = cellfun(@(x) x(1:end-4), ms1_names, 'UniformOutput', false);

% Suffixes to ignore
suffixes = {'_HCDFT', '_HCDIT', '_CIDFT', '_CIDIT'};

% Check whether each .mgf file has a corresponding .ms1 file
for i = 1:length(mgf_names)
    current_mgf_base = mgf_names{i};
    
    % Try to strip known suffixes strictly from the end
    cleaned_name = current_mgf_base;
    for s = 1:length(suffixes)
        suffix = suffixes{s};
        if endsWith(cleaned_name, suffix)
            cleaned_name = cleaned_name(1:end - length(suffix));
            break; % Assume only one suffix applies
        end
    end
    
    % Find matches in ms1 files
    % 1. Try to match with the cleaned name (suffix removed)
    match_cleaned = strcmp(cleaned_name, ms1_names);
    
    % 2. If cleaned name is different from original (i.e. suffix was removed),
    %    also try to match with the original name (full match)
    match_original = false(size(ms1_names));
    if ~strcmp(cleaned_name, current_mgf_base)
        match_original = strcmp(current_mgf_base, ms1_names);
    end

    % Combine matches: valid if EITHER cleaned name OR original name matches
    matches = match_cleaned | match_original;
    num_matches = sum(matches);
    
    if num_matches == 0
        error('File %s.mgf does not have a corresponding .ms1 file (expected base: %s or %s) in %s.', mgf_names{i}, cleaned_name, current_mgf_base, obj.m_specPath);
    elseif num_matches > 1
        warning('File %s.mgf corresponds to multiple .ms1 files. Result might be ambiguous.', mgf_names{i});
    end
end
end