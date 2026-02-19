function [obj, created] = ensureMgfDatasetIO(obj)
% Ensure MGF dataset IO is initialized
% Output:
%   obj (CMSMSPepDeconv)
%       updated instance
%   created (1 x 1 logical)
%       true if created in this call

created = false;
if isempty(obj.m_cMgfDatasetIO)
    obj.m_cMgfDatasetIO = CMgfDatasetIO;
    obj.m_cMgfDatasetIO.Init(obj.m_specPath);
    obj.m_cMgfDatasetIO.SetMap();
    obj.m_cMgfDatasetIO.SetFidmap();
    created = true;
end
end
