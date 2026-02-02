function build_mapping(obj)
    % Build mgf->ms1 mapping based on file stems and suffix rules
    mgf_files = dir(fullfile(obj.m_specPath, '*.mgf'));
    ms1_files = dir(fullfile(obj.m_specPath, '*.ms1'));

    if isempty(mgf_files)
        % It is acceptable to have no files during some initializations, 
        % but if we are mapping, it usually implies we expect files.
        % For now, just return.
        return;
    end

    % Get stems
    mgf_names = {mgf_files.name};
    ms1_names = {ms1_files.name};
    
    % Remove extension to get "stem"
    mgf_stems = cellfun(@(x) x(1:end-4), mgf_names, 'UniformOutput', false);
    ms1_stems = cellfun(@(x) x(1:end-4), ms1_names, 'UniformOutput', false);

    for i = 1:length(mgf_stems)
        current_mgf_base = mgf_stems{i};
        
        % Try to strip known suffixes strictly from the end
        cleaned_name = current_mgf_base;
        for s = 1:length(obj.m_suffixes)
            suffix = obj.m_suffixes{s};
            if endsWith(cleaned_name, suffix)
                cleaned_name = cleaned_name(1:end - length(suffix));
                break; % Only one suffix
            end
        end
        
        % Find matches in ms1 files
        % 1. Try to match with the cleaned name (suffix removed)
        match_cleaned = strcmp(cleaned_name, ms1_stems);
        
        % 2. If cleaned name is different from original (i.e. suffix was removed),
        %    also try to match with the original name (full match)
        match_original = false(size(ms1_stems));
        if ~strcmp(cleaned_name, current_mgf_base)
            match_original = strcmp(current_mgf_base, ms1_stems);
        end

        % Combine matches: valid if EITHER cleaned name OR original name matches
        matches = match_cleaned | match_original;
        num_matches = sum(matches);
        idxs = find(matches);

        if num_matches == 0
            if strcmp(cleaned_name, current_mgf_base)
                 error('File %s.mgf (expected base: %s) does not have a corresponding .ms1 file in %s.', current_mgf_base, cleaned_name, obj.m_specPath);
            else
                 error('File %s.mgf (expected base: %s or %s) does not have a corresponding .ms1 file in %s.', current_mgf_base, cleaned_name, current_mgf_base, obj.m_specPath);
            end
        elseif num_matches > 1
            warning('File %s.mgf corresponds to multiple .ms1 files. Result might be ambiguous. Using %s.', current_mgf_base, ms1_stems{idxs(1)});
        end

        if num_matches >= 1
            obj.m_mgf2ms1_map(current_mgf_base) = ms1_stems{idxs(1)};
        end
        
    end
end