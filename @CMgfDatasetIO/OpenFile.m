function fid=OpenFile(~,strFilename)
% Open file with error checking, returns an error if opening fails

fid=fopen(strFilename,'r');
if 0>=fid
    error('Filed to open dataset file!');
end

end
