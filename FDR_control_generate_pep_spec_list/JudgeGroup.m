function [DecoyType,GroupType,OriginalPeptide,scores,numrst,I] = JudgeGroup(result,TagType,DecoyTag,GroupTag)
% input:
%           result (1 x N struct)
%               result read from the search engine result file
%           TagType (1 x 1 char/string)
%               grouping type (Protein/Modification/SNP)
%           DecoyTag (1 x 1 char/string)
%               tag indicating decoy match
%           GroupTag (char/string or cell)
%               keyword(s) for grouping
%           databaseHash
%               only used in SNP problem
% output:
%           DecoyType (1 x N double/logical)
%               decoy indicator per PSM
%           GroupType (1 x N double/logical)
%               out-of-group indicator per PSM (0 is wanted)
%           OriginalPeptide (N x 1 cell)
%               only used in SNP problem
%           scores (1 x N double)
%               search engine scores (sorted, descending)
%           numrst (1 x 1 double/int)
%               number of results
%           I (1 x N double/int)
%               indices of PSMs in descending sort

Scores = [result(:).Score];
[scores,I] = sort(Scores,'descend');
numrst = numel(I);% Calculate the number of elements in I
Result = result(I);
DecoyType = zeros(1,length(Result));
GroupType = zeros(1,length(Result));
OriginalPeptide = cell(length(Result),1);

for i = 1:length(Result)
    %%%%%%%%%%%%%%%%%%%%%%%%%%Judge Target and Decoy type%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if contains(upper(Result(i).protein),upper(DecoyTag))
        DecoyType(i) = 0;
    else
        DecoyType(i) = 1;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%END%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%Judge Tag type: Protein or Modification%%%%%%%%%%%%%%%%%%%%
    if strcmpi(TagType,'Protein')
        if size(GroupTag,1) == 1 % only one tag
            if strfind(Result(i).protein,GroupTag)
                GroupType(i) = 0;
            else
                GroupType(i) = 1;
            end
        else % many tags, in string cells format
            for idx_gt = 1:size(GroupTag,1)
                if contains(Result(i).protein,GroupTag{idx_gt})
                    GroupType(i) = 0;
                    break;
                else
                    GroupType(i) = 1;
                end
            end
        end
    else % strcmpi(TagType,'Modification')
        if strfind(Result(i).modification,GroupTag)
            GroupType(i) = 0;
        else
            GroupType(i) = 1;
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%END%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

end



