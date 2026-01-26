function check_whether_ms12_mgf_name_match(obj)
% Check whether the name of ms1/ms2 is matched with mgf
% DEPRECATED: This logic has been moved to CMsFileMapper class.
% Keeping this method as a wrapper for backward compatibility or direct calls.
% TODO: Remove this method if no longer needed.

if isempty(obj.m_cMsFileMapper)
    obj.m_cMsFileMapper = CMsFileMapper(obj.m_specPath);
else
    % If it exists, force check by rebuilding mapping or assume it's done.
    obj.m_cMsFileMapper.build_mapping();
end

end