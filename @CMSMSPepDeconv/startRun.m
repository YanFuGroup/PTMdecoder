function obj = startRun(obj, is_record_fragment_information)
% Start the run for MSMS level and Peptide level IMP discrimination and quantification
% Input:
%   obj (CMSMSPepDeconv)
%       Processor instance
%   is_record_fragment_information (1 x 1 logical, optional)
%       Whether to record fragment information
% Output:
%   obj (CMSMSPepDeconv)
%       Updated instance

if nargin < 2
    is_record_fragment_information = false;
end

% Initialize shared resources lazily
[obj, ~] = obj.ensureMsFileMapper();
[obj, mgf_created_here] = obj.ensureMgfDatasetIO();
if mgf_created_here
    cleanup_mgf = onCleanup(@() obj.m_cMgfDatasetIO.CloseAllFile()); %#ok<NASGU>
end

% Quantification each IMP for each PSM
obj = obj.MSMSLevelRun(is_record_fragment_information);

% Quantification each modified peptide according to each PSMs
obj = obj.quant_IMP();

end