function [cStr]=massArraTostr(obj,massArrangement,fixedPosMod,dictVariMod,inxSites)
% Organize modification combinations in string form
% Input:
%   massArrangement (M x S double) - the various mass arrangements of modifications on the peptide
%   fixedPosMod (K x 3 cell) - the fixed modification [modification type-specificity-modification mass] matrix
%   dictVariMod (R x 3 cell) - the variable modification [modification type-specificity-modification mass] matrix, the order is consistent with all modification types specified by the user
%   inxSites (1 x S double/int) - the amino acid positions where modifications occur on the peptide, arranged in ascending order
% Output: 
%   cStr (M x 1 cell) - IMP strings, which is a column cell array,
%       in the form of _{Propionyl}GK{Acetyl}GGK{Propionyl}GLGK{Propionyl}GGAK{Propionyl}R_
cStr=cell(size(massArrangement,1),1);
for t=1:size(massArrangement,1)% by row
    strExtPep=['_',obj.m_pepSeq,'_'];
    strExtPep = cellstr(strExtPep');
    % fixed modification
    for idxFixMod = 1:size(fixedPosMod,1)
        strExtPep{fixedPosMod{idxFixMod,1}+1} = [strExtPep{fixedPosMod{idxFixMod,1}+1},...
            '{',fixedPosMod{idxFixMod,2},'}'];
    end
    % variable modification
    for s=1:length(inxSites)% Iterate through sites
        mass=massArrangement(t,s);
        for r=1:length(dictVariMod)
            if dictVariMod{r,3}==mass% Find the name of the modification by mass
                strExtPep{inxSites(s)+1} = [strExtPep{inxSites(s)+1},'{',dictVariMod{r,1},'}'];
                break;
            end
        end
    end
    cStr{t}=[strExtPep{:}];
end
end