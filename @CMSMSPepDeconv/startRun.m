function obj = startRun(obj, is_record_fragment_information)
% Start the run for MSMS level and Peptide level IMP discrimination and  quantification

if nargin < 2
    is_record_fragment_information = false;
end

% Check whether the name of ms1/ms2 is matched with mgf
% obj.check_whether_ms12_mgf_name_match();
% Initialize the file mapper
if isempty(obj.m_cMsFileMapper)
    obj.m_cMsFileMapper = CMsFileMapper(obj.m_specPath);
end

% Indexing the mgf
obj.m_cMgfDatasetIO = CMgfDatasetIO;
obj.m_cMgfDatasetIO.Init(obj.m_specPath);
obj.m_cMgfDatasetIO.SetMap();
obj.m_cMgfDatasetIO.SetFidmap();

% Quantification each IMP for each PSM
obj = obj.MSMSLevelRun(is_record_fragment_information);

% Quantification each modified peptide according to each PSMs
obj = obj.PepLevelRun();

obj.m_cMgfDatasetIO.CloseAllFile();

end