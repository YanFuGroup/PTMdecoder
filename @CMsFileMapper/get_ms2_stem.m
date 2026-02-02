function ms2_stem = get_ms2_stem(obj, mgf_stem)
    % Returns the stem of the MS2 file corresponding to the given MGF stem
    % Input:
    %   mgf_stem (1 x 1 char/string)
    %       mgf file stem (without extension)
    % Output:
    %   ms2_stem (1 x 1 char/string)
    %       ms2 file stem
    ms1_stem = obj.get_ms1_stem(mgf_stem);
    ms2_stem = strrep(ms1_stem, '.ms1', '.ms2');
end