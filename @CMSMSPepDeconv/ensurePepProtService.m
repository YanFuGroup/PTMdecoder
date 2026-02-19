function [obj, created] = ensurePepProtService(obj)
% Ensure peptide-protein service is initialized
% Output:
%   obj (CMSMSPepDeconv)
%       updated instance
%   created (1 x 1 logical)
%       true if created in this call

created = false;
if isempty(obj.CPepProtService)
    obj.CPepProtService = CPepProtService(obj.m_fastaFile, obj.m_regular_express, obj.m_filtered_res_file_path);
    created = true;
end
end
