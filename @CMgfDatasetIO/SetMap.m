function SetMap(obj)
% Create an index from spectrum name to position in mgf files for all spectra.
% Input:
%   obj (CMgfDatasetIO)
%       dataset IO instance

dataset_files = dir(fullfile(obj.m_strFoldname,'*.mgf'));
dataset_names = {dataset_files.name}'; % Arrange file names in a column

fprintf('%s%d\n','Total Number of Dataset:',length(dataset_files));

% begin building up
for i = 1:length(dataset_names) % Iterate over each spectral dataset
    filename = fullfile(obj.m_strFoldname,dataset_names{i,1});
    fprintf('%s','Spectra: ');
    nSpec = 0;
    
    % load mgf map when existing a intermediate map file
    [mgfpath,mgfname] = fileparts(filename);
    mgf_mapfile = fullfile(mgfpath,[mgfname,'_MGF_map.mat']);
    ct_prt = 0; % length of print string
    if exist(mgf_mapfile,"file")
        load(mgf_mapfile,'mgf_map')
        nSpec = mgf_map.Count;
    else
        fid = fopen(filename,'r');
        if 0>=fid
            disp(['Failed to open dataset: ',filename,'!']);
            return
        end

        mgf_map = containers.Map();

        while ~feof(fid)
            strLine = fgets(fid);

            % the begin of one spectrum
            if strncmp(strLine,'BEGIN IONS',10)
                nSpec = nSpec+1;

                if mod(nSpec,2000) == 0
                    fprintf(repmat('\b',1,ct_prt));
                    ct_prt = fprintf('%d',nSpec);
                end

                % The starting position of the file corresponding to each spectrum is BEGIN IONS
                iPosition = ftell(fid)-length(strLine);

                while(~strncmp(strLine,'TITLE=',6)) % go to TITLE
                    strLine = fgetl(fid);
                end

                strSpecName = strLine(7:end);
                mgf_map(strSpecName) = iPosition;

                if contains(strSpecName,'.')
                    strScanNum = regexp(strSpecName,'\.','split');
                    MS2ScanI = str2double(strScanNum{2});
                    if ~isnan(MS2ScanI)
                        mgf_map(strScanNum{2}) = iPosition;
                    end
                end
            end
        end
        save(mgf_mapfile,'mgf_map');
        fclose(fid);
    end
    
    % save spectral index for each dataset
    obj.m_mapDatasetIdx(dataset_names{i,1}) = mgf_map;
    fprintf(repmat('\b',1,ct_prt));
    fprintf('%d\n',nSpec);
end

end