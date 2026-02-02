function ms1_stem = get_ms1_stem(obj, mgf_stem)
    % Returns the stem of the MS1 file corresponding to the given MGF stem
    % Input:
    %   mgf_stem (1 x 1 char/string)
    %       mgf file stem (without extension)
    % Output:
    %   ms1_stem (1 x 1 char/string)
    %       ms1 file stem
    if isKey(obj.m_mgf2ms1_map, mgf_stem)
        ms1_stem = obj.m_mgf2ms1_map(mgf_stem);
    else
        % If not found (unlikely if build_mapping ran and errored on missing), 
        % maybe the file was added later or the check logic was bypassed?
        % Or maybe the user queried a name that simply does not exist.
        error('MGF stem "%s" not found in mapping.', mgf_stem);
    end
end