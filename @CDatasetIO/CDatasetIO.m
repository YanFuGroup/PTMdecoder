classdef CDatasetIO<handle
% Dataset IO class
    properties(SetAccess=protected)
        m_strFoldname;  % Directory name
    end
    
    methods (Access=protected)
        SetMap(obj);
        % Input:
        %   obj (CDatasetIO)
        %       dataset IO base instance
    end
end