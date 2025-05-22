function CloseAllFile(obj)
% Close all file in the mgf dictionary
Fids=values(obj.m_mapFid);
for i=1:length(Fids)
    fclose(Fids{i});
end

end