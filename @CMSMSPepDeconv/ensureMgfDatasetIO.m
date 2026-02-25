function obj = ensureMgfDatasetIO(obj)
% Ensure MGF dataset IO is initialized
% Output:
%   obj (CMSMSPepDeconv)
%       updated instance

if isempty(obj.m_cMgfDatasetIO)
    obj.m_cMgfDatasetIO = CMgfDatasetIO(obj.m_specPath);
end
end
