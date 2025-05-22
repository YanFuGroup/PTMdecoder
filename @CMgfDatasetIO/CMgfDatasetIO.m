classdef CMgfDatasetIO<CDatasetIO
% MGF spectrum data IO
    properties
        m_mapFid;       % Dictionary: from MGF filename to fid
        m_mapDatasetIdx;% Dictionary: from MGF name to a dictionary mapping spectrum title to its position in mgf file
    end
    methods  
        Init(obj,cMLocParam); % Initialization
        
        SetFidmap(obj,strFoldname); % Load m_mapFid from directory strFoldname, which contains MGF files
        
        fid=OpenFile(obj,strFilename); % Given a filename, returns a file pointer; prompts if opening fails
        
        CloseFile(obj,fid); % Close the specified fid
        
        CloseAllFile(obj); % Close all files in the MGF dictionary
        
        SetMap(obj); % Load m_mapDatasetIdx from the MGF files in the directory 

        % read the specified spectrum in the mgf file
        [Peaks,Charge,PrecursorMZ]=read_oneSpec(obj,filename,specname);
    end
end
