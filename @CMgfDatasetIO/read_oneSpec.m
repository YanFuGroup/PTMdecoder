function [Peaks,Charge,PrecursorMZ]=read_oneSpec(obj,filename,specname)
% Read one spectrum from the dataset and output peak information
% Input:
%   obj (CMgfDatasetIO)
%       dataset IO instance
%   filename (1 x 1 char/string)
%       dataset name (mgf file name)
%   specname (1 x 1 char/string)
%       spectrum name
% Output:
%   Peaks (N x 2 double)
%       peak list [m/z, intensity]
%   Charge (1 x 1 double/int)
%       charge state
%   PrecursorMZ (1 x 1 double)
%       precursor ion m/z
% ATTENTION: Before using this function, call SetFidmap so files can be opened.

BUF_LENGTH=500;
Peaks=zeros(BUF_LENGTH,2);% Pre-allocate
nPeaks=0;
mapSpecIdx=obj.m_mapDatasetIdx(filename);% mapSpecIdx is a dictionary from spectrum name (title) to its position in the file
iPosition=mapSpecIdx(specname);% File position
fid=obj.m_mapFid(filename);
if -1==fid
    error(['Failed to open file ' filename]);
end
status=fseek(fid,iPosition,'bof');% Move fid to that position, and check the 
if 0~=status
    error('Failed to locate spectral position!');
end

% Read information from the spectrum
strLine=fgetl(fid);% Skip the first line
while~strcmp(strLine,'END IONS')% Keep getting peaks as long as it's not END IONS
    % peaks
    strLine=fgetl(fid);
    if isstrprop(strLine(1),'digit')% If the first character is a digit, count it and record the peak m/z and intensity
        nPeaks=nPeaks+1;
        Peaks(nPeaks,:)=sscanf(strLine,'%f\t%f',[1,2]);
        continue;
    else
        if startsWith(strLine,'CHARGE=')
            Charge=sscanf(strLine,'CHARGE=%d+');% Only supports positive charge mode
        elseif startsWith(strLine,'PEPMASS=')
            PrecursorMZ=sscanf(strLine,'PEPMASS=%f');
        end
    end
end
Peaks(nPeaks+1:end,:)=[];
% fclose(fid);
end