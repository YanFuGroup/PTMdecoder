classdef CDatasetIO<handle
% Dataset IO class
    properties(SetAccess=protected)
        m_strFoldname;  % Directory name
    end
    
    methods
        Init(obj,cMLocParam);
        
        SetMap(obj);
    end
end