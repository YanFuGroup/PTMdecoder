function [result] = ReadDatResultFolder(path)

result = [];
datafiles = dir(fullfile(path,'*.dat'));
nfiles = length(datafiles);
for i=1:nfiles
    plfilepath = fullfile(path,datafiles(i).name);
    fprintf('Reading %s...',plfilepath);
    [resulti] = ReadDatResult(plfilepath);   
    result = [result resulti];
    fprintf('\n')   
end

end