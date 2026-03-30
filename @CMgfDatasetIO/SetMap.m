function SetMap(obj)
% Create an index from spectrum name to position in mgf files for all spectra.
% Input:
%   obj (CMgfDatasetIO)
%       dataset IO instance

dataset_files = dir(fullfile(obj.m_strFoldname,'*.mgf'));
dataset_names = {dataset_files.name}'; % Arrange file names in a column

CLogger.info('Total number of datasets: %d', length(dataset_files));

% begin building up
for i = 1:length(dataset_names) % Iterate over each spectral dataset
    filename = fullfile(obj.m_strFoldname,dataset_names{i,1});
    dataset_size = dataset_files(i).bytes;
    progress_label = sprintf('MGF map %s', dataset_names{i,1});
    nSpec = 0;
    
    % load mgf map when existing a intermediate map file
    [mgfpath,mgfname] = fileparts(filename);
    mgf_mapfile = fullfile(mgfpath,[mgfname,'_MGF_map.mat']);
    if exist(mgf_mapfile,"file")
        load(mgf_mapfile,'mgf_map')
        nSpec = mgf_map.Count;
        CLogger.debug('Loaded existing MGF map: %s', mgf_mapfile);
    else
        fid = fopen(filename,'r');
        if 0>=fid
            CLogger.warn('Failed to open dataset: %s', filename);
            return
        end

        mgf_map = containers.Map();

        while ~feof(fid)
            strLine = fgets(fid);

            % the begin of one spectrum
            if strncmp(strLine,'BEGIN IONS',10)
                nSpec = nSpec+1;

                if mod(nSpec,2000) == 0
                    CLogger.progress(progress_label, ftell(fid), max(dataset_size, 1));
                end

                % The starting position of the file corresponding to each spectrum is BEGIN IONS
                iPosition = ftell(fid)-length(strLine);

                while(~strncmp(strLine,'TITLE=',6)) % go to TITLE
                    strLine = fgetl(fid);
                end

                strSpecName = strLine(7:end);
                mgf_map(strSpecName) = iPosition;

                try
                    MS2ScanI = CMS2SpecNameUtils.parseMS2ScanNumber(strSpecName);
                    scan_key = num2str(MS2ScanI);
                    if ~isequal(scan_key, strSpecName)
                        mgf_map(scan_key) = iPosition;
                    end
                catch
                    % Keep raw TITLE mapping only when scan number is not parseable.
                end
            end
        end
        CLogger.progress(progress_label, max(dataset_size, 1), max(dataset_size, 1));
        save(mgf_mapfile,'mgf_map');
        fclose(fid);
    end
    
    % save spectral index for each dataset
    obj.m_mapDatasetIdx(dataset_names{i,1}) = mgf_map;
    CLogger.info('Indexed dataset %s with %d spectra.', dataset_names{i,1}, nSpec);
end

end