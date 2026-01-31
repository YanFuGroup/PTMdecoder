function fid=OpenFile(~,strFilename)
% Open file with error checking, returns an error if opening fails
% Input:
%   strFilename (1 x 1 char/string)
%       file path
% Output:
%   fid (1 x 1 double/int)
%       file identifier

fid=fopen(strFilename,'r');
if 0>=fid
    error('Filed to open dataset file!');
end

end
