classdef CMgfDatasetIO<CDatasetIO
% MGF spectrum data IO
    properties
        m_mapFid;       % Dictionary: from MGF filename to fid
        m_mapDatasetIdx;% Dictionary: from MGF name to a dictionary mapping spectrum title to its position in mgf file
    end
    methods
        function obj = CMgfDatasetIO(strDatasetFoldname)
        % Construct MGF dataset IO instance.
        % Input:
        %   strDatasetFoldname (1 x 1 char/string, optional)
        %       dataset folder path; if provided, initialization is completed directly.
            if nargin < 1
                return;
            end
            if isempty(strDatasetFoldname)
                error('strDatasetFoldname is required.');
            end

            % Close old file handles when re-initializing the same object
            if ~isempty(obj.m_mapFid)
                obj.CloseAllFile();
            end

            obj.m_mapDatasetIdx=containers.Map();% Create an instance of containers.Map
            obj.m_mapFid=containers.Map();% Create an instance of containers.Map
            obj.m_strFoldname=strDatasetFoldname;

            % One-step setup: build spectrum index and open file handles
            obj.SetMap();
            obj.SetFidmap();
        end

        delete(obj); % Destructor: close all opened MGF files
        
        fid=OpenFile(obj,strFilename); % Given a filename, returns a file pointer; prompts if opening fails
        
        CloseFile(obj,fid); % Close the specified fid
        
        CloseAllFile(obj); % Close all files in the MGF dictionary

        % read the specified spectrum in the mgf file
        [Peaks,Charge,PrecursorMZ]=read_oneSpec(obj,filename,specname);
    end

    methods (Access=protected)
        SetMap(obj); % Load m_mapDatasetIdx from the MGF files in the directory
    end

    methods (Access=private)
        SetFidmap(obj); % Load m_mapFid from obj.m_strFoldname, which contains MGF files
    end
end
