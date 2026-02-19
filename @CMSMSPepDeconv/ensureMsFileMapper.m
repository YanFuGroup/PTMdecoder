function [obj, created] = ensureMsFileMapper(obj)
% Ensure MGF<->MS1 file mapper is initialized
% Output:
%   obj (CMSMSPepDeconv)
%       updated instance
%   created (1 x 1 logical)
%       true if created in this call

created = false;
if isempty(obj.m_cMsFileMapper)
    obj.m_cMsFileMapper = CMsFileMapper(obj.m_specPath);
    created = true;
end
end
