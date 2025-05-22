function obj = Init(obj)
% Initialize, read sequence from file to m_strOneLine

fin = fopen(obj.m_strFilePath,'r');
if 0>=fin
    error(['Can not open file: ',obj.m_strFilePath]);
end

rows = 0;
while ~feof(fin)
    rows = rows+sum(fread(fin,10000,'char')==10);
end
cstrFileContents = repmat({''},1,ceil(rows/2)+2);
cstrProtName = repmat({''},1,ceil(rows/2)+2);

idx_content = 0;
idx_prot_name = 0;
frewind(fin);
while ~feof(fin)
    strline = fgetl(fin);
    if isempty(strline)
        continue;
    end

    idx_content = idx_content+1;
    % If it starts with ">", it is the beginning of a protein, replace this line with $
    if strline(1)=='>'
        cstrFileContents{idx_content} = '$';
        match = regexp(strline, obj.m_regular_exp, 'tokens', 'once');
        idx_prot_name = idx_prot_name + 1;
        cstrProtName{idx_prot_name} = match{1};
    else
        cstrFileContents{idx_content} = strline;
    end
end
cstrFileContents{idx_content+1} = '$';
cstrFileContents(idx_content+2:end) = [];
cstrProtName(idx_prot_name+1:end) = [];

obj.m_strOneLine = cell2mat(cstrFileContents);
obj.m_strProtName = cstrProtName;
obj.m_strProtPos = strfind(obj.m_strOneLine,'$') + 1;
obj.m_strProtPos(end) = [];
fclose(fin);
end
