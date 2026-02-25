function CloseAllFile(obj)
% Close all file in the mgf dictionary
% Input:
%   obj (CMgfDatasetIO)
%       dataset IO instance
if isempty(obj.m_mapFid)
    return;
end

Fids=values(obj.m_mapFid);
for i=1:length(Fids)
    fid = Fids{i};
    if isnumeric(fid) && isscalar(fid) && fid > 0
        fclose(fid);
    end
end

end