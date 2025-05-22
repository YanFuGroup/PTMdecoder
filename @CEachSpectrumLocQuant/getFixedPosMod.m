function [fixedPosMod]=getFixedPosMod(obj)
% According to the user-specified fixed modification, find their positions on the peptide sequence and record the information
% Output: 
%   fixedPosMod
%       if none then empty, 
%       if present each row is a fixed modification, each row is [position, name, mass], such as: [2, 'Acetyl[K]', 42.0106]
if isempty(obj.m_fixedModNameMass)
    fixedPosMod=[];
    return;
end
fixedPosMod=cell(length(obj.m_pepSeq)+2,3);
% Add fixed modifications
for idx=1:size(obj.m_fixedModNameMass,1)
    strSpecfct=obj.m_fixedModNameMass{idx,2};   % string of specificity

    % Only consider the case of a single amino acid (can be N/C terminus), not multiple amino acids in parallel
    if contains(strSpecfct,'N-term')
        cellTemp=split(strSpecfct,'N-term');
        if isequal(cellTemp{1},'Protein') && ~obj.m_isProtN
            continue;
        end
        if ~isempty(cellTemp{2})
            if ~startsWith(obj.m_pepSeq,cellTemp{2})
                continue; % N-terminal specificity, with amino acid restriction
            else
                fixedPosMod{1,1} = 1;
                fixedPosMod{1,2} = obj.m_fixedModNameMass{idx,1};
                fixedPosMod{1,3} = obj.m_fixedModNameMass{idx,3};
            end
        else % N-terminus without specified amino acid
            fixedPosMod{1,1} = 0;
            fixedPosMod{1,2} = obj.m_fixedModNameMass{idx,1};
            fixedPosMod{1,3} = obj.m_fixedModNameMass{idx,3};
        end
    elseif contains(strSpecfct,'C-term')
        cellTemp=split(strSpecfct,'C-term');
        if isequal(cellTemp{1},'Protein') && ~obj.m_isProtC
            continue;
        end
        if ~isempty(cellTemp{2})
            if ~endsWith(obj.m_pepSeq,cellTemp{2})
                continue; % C-terminal specificity, with amino acid restriction
            else
                fixedPosMod{length(obj.m_pepSeq)+2,1} = length(obj.m_pepSeq);
                fixedPosMod{length(obj.m_pepSeq)+2,2} = obj.m_fixedModNameMass{idx,1};
                fixedPosMod{length(obj.m_pepSeq)+2,3} = obj.m_fixedModNameMass{idx,3};
            end
        else % C-terminus without specified amino acid
            fixedPosMod{length(obj.m_pepSeq)+2,1} = length(obj.m_pepSeq)+1;
            fixedPosMod{length(obj.m_pepSeq)+2,2} = obj.m_fixedModNameMass{idx,1};
            fixedPosMod{length(obj.m_pepSeq)+2,3} = obj.m_fixedModNameMass{idx,3};
        end
    else
        % Only one amino acid
        vecPos=strfind(obj.m_pepSeq,strSpecfct);
        for iPos=1:length(vecPos)
            fixedPosMod{vecPos(iPos)+1,1} = vecPos(iPos);
            fixedPosMod{vecPos(iPos)+1,2} = obj.m_fixedModNameMass{idx,1};
            fixedPosMod{vecPos(iPos)+1,3} = obj.m_fixedModNameMass{idx,3};
        end
    end
end
fixedPosMod(cellfun(@isempty,fixedPosMod(:,1)),:)=[];
end