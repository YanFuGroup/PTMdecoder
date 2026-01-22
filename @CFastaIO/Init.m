function obj = Init(obj)
% Initialize, read sequence from file to m_mapProt

% fin = fopen(obj.m_strFilePath,'r');
% if 0>=fin
%     error(['Can not open file: ',obj.m_strFilePath]);
% end
% 
% rows = 0;
% while ~feof(fin)
%     rows = rows+sum(fread(fin,100000,'char')==10);
% end
% cstrFileContents = repmat({''},1,ceil(rows/2)+2);
% cstrProtName = repmat({''},1,ceil(rows/2)+2);
% 
% idx_content = 0;
% idx_prot_name = 0;
% frewind(fin);
% while ~feof(fin)
%     strline = fgetl(fin);
%     if isempty(strline)
%         continue;
%     end
% 
%     idx_content = idx_content+1;
%     % If it starts with ">", it is the beginning of a protein, replace this line with $
%     if strline(1)=='>'
%         cstrFileContents{idx_content} = '$';
%         match = regexp(strline, obj.m_regular_exp, 'tokens', 'once');
%         idx_prot_name = idx_prot_name + 1;
%         cstrProtName{idx_prot_name} = match{1};
%     else
%         cstrFileContents{idx_content} = strline;
%     end
% end
% cstrFileContents{idx_content+1} = '$';
% cstrFileContents(idx_content+2:end) = [];
% cstrProtName(idx_prot_name+1:end) = [];
% 
% obj.m_strOneLine = cell2mat(cstrFileContents);
% obj.m_strProtName = cstrProtName;
% obj.m_strProtPos = strfind(obj.m_strOneLine,'$') + 1;
% obj.m_strProtPos(end) = [];
% fclose(fin);

% New implementation: use containers.Map to store protein name and sequence
fin = fopen(obj.m_strFilePath, 'r');
if fin <= 0
    error(['Can not open file: ', obj.m_strFilePath]);
end

obj.m_mapProt = containers.Map('KeyType', 'char', 'ValueType', 'char');

currentProtName = '';
currentSeq = '';

while ~feof(fin)
    strLine = strtrim(fgetl(fin));
    if isempty(strLine)
        continue;
    end
    
    if strLine(1) == '>'
        % Save previous protein if it exists
        if ~isempty(currentProtName)
            obj.m_mapProt(currentProtName) = currentSeq;
        end
        
        % Extract protein name using regular expression
        match = regexp(strLine, obj.m_regular_exp, 'tokens', 'once');
        if ~isempty(match)
            currentProtName = match{1};
        else
            warning(['Protein name not found with regular expression in header: ', strLine]);
            currentProtName = strLine(2:end); % Fallback to rest of header
        end
        currentSeq = '';
    else
        % Concatenate sequence lines
        currentSeq = [currentSeq, strLine];
    end
end

% Save the last protein
if ~isempty(currentProtName)
    obj.m_mapProt(currentProtName) = currentSeq;
end

fclose(fin);
end
