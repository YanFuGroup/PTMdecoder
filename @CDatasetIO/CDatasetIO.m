classdef CDatasetIO<handle
% Dataset IO class
    properties(SetAccess=protected)
        m_strFoldname;  % Directory name
    end
    
    methods
        Init(obj,cMLocParam);
        % Input:
        %   obj (CDatasetIO)
        %       dataset IO base instance
        %   cMLocParam (1 x 1 char/string)
        %       dataset folder path
        
        SetMap(obj);
        % Input:
        %   obj (CDatasetIO)
        %       dataset IO base instance
    end
end