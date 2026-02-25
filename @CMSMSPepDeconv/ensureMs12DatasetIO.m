function [obj, created] = ensureMs12DatasetIO(obj)
% Ensure MS1/MS2 dataset IO is initialized
% Output:
%   obj (CMSMSPepDeconv)
%       updated instance
%   created (1 x 1 logical)
%       true if created in this call

created = false;
if isempty(obj.m_cMs12DatasetIO)
    obj.m_cMs12DatasetIO = CMS12DatasetIO(obj.m_specPath, obj.m_ms1_tolerance);
    created = true;
end
end
