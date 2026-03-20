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
% ATTENTION: Before using this function, call Init so index and file handles are ready.

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

% Initialize
Charge = 0;
PrecursorMZ = 0;

% Read the spectrum
fgetl(fid);% Skip the first line
while ~feof(fid)
    strLine=fgetl(fid);

    if isempty(strLine)
        continue;
    end

    % End of the spectrum ahead, return empty peaks in advance
    if strcmp(strLine, 'END IONS')
        Peaks = zeros(0, 2);
        break;
    end

    firstChar = strLine(1);
    if (firstChar >= '0' && firstChar <= '9') || firstChar == '.'
        firstPeak = sscanf(strLine, '%f', 2);   % Parse the 'this' peak line
        peakData = textscan(fid, '%f %f');      % Parse the following peak lines until a non-peak line is encountered
        Peaks = [firstPeak'; peakData{1}, peakData{2}];
        break;
    else
        if startsWith(strLine,'CHARGE=')
            Charge=sscanf(strLine,'CHARGE=%d+');% Only supports positive charge mode
        elseif startsWith(strLine,'PEPMASS=')
            PrecursorMZ=sscanf(strLine,'PEPMASS=%f');
        end
    end
end
end