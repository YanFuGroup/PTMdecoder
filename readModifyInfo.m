% Create a dictionary mapping modification names to mono masses using the pFind studio style modification library file
% Input: strPathname (1 x N char/string) - the file name
% Output: mapModification (containers.Map) - a dictionary mapping modification names (with amino acid specificity) to mono mass offsets
function mapModification=readModifyInfo(strPathname)
% read modification informations
mapModification=containers.Map();

fid=fopen(strPathname,'r');
if 0>=fid
    error('Failed to open modify.ini file!');
end

% Skip the first line, keep the part after the equal sign on the second line, and convert it to a numeric type
fgetl(fid);
strLine=fgetl(fid);
[~,strEnd]=strtok(strLine,'=');
nModifications=str2double(strEnd(2:end));% Number of modifications read from the second line

nModify=0;
while ~feof(fid)
    strLine=fgetl(fid);
    
    if isempty(strLine)
        continue;
    end
    
    if strncmpi(strLine,'name',4)% If it is a 'name=' line, then read the next line, which is its detailed information line
        strLine=fgetl(fid);
        
        [strName,strValue]=strtok(strLine,'=');
        
        Idx=strfind(strValue,' ');
        strValue=str2double(strValue(Idx(2)+1:Idx(3)-1));% The mono mass of the modification is between the second and third blank space ' ', convert it to a numeric type
        
        mapModification(strName)=strValue;% Create a dictionary mapping modification names to mono masses
        
        nModify=nModify+1;% Record the number of modifications in the dictionary
    end
end

if nModifications~=nModify
    % Whether the number of modifications read is consistent with the one declared at the beginning of the file
    error('Failed to read modify.ini file! The number of read modification is not consistent with the header.');
end

mapModification('null') = 0;
mapModification('NULL') = 0;
mapModification('Null') = 0;
mapModification('unknown') = 0;
fclose(fid);
end

