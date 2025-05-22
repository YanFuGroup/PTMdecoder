function [DecoyType,GroupType,OriginalPeptide,scores,numrst,I] = JudgeGroup(result,TagType,DecoyTag,GroupTag)
% input:
%           result
%               the result read from the search engine result file
%           TagType
%               the tag showing the type of grouping problem,
%               Protein, Modification, SNP and so on
%           DecoyTag
%               the tag indicating the decoy match
%           GroupTag
%               the key word (tag) of grouping
%           databaseHash
%               only used in SNP problem
% output:
%           DecoyType
%               one for each PSM, indicating whether it is a decoy match
%           GroupType
%               one for each PSM, indicating whether it is out of the
%               target group, 0 is wanted
%           OriginalPeptide
%               only used in SNP problem
%           scores
%               the score of search engine, sorted
%           numrst
%               the number of results
%           I
%               the index of the PSMs in descend sort

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



