function [fixedPosMod]=getFixedPosMod(ctx)
% According to the user-specified fixed modification, find their positions on the peptide sequence and record the information
% Input:
%   ctx (struct)
%       Required fields: m_pepSeq, m_fixedModNameMass, m_isProtN, m_isProtC.
% Output:
%   fixedPosMod (K x 3 cell)
%       Each row: [position, mod_name, mod_mass], such as: [2, 'Acetyl[K]', 42.0106]. Empty if no fixed mods.

if isempty(ctx.m_fixedModNameMass)
	fixedPosMod=[];
	return;
end
fixedPosMod=cell(length(ctx.m_pepSeq)+2,3);
for idx=1:size(ctx.m_fixedModNameMass,1)
	strSpecfct=ctx.m_fixedModNameMass{idx,2};

	% Only consider the case of a single amino acid (can be N/C terminus), not multiple amino acids in parallel
	if contains(strSpecfct,'N-term')
		cellTemp=split(strSpecfct,'N-term');
		if isequal(cellTemp{1},'Protein') && ~ctx.m_isProtN
			continue;
		end
		if ~isempty(cellTemp{2})
			if ~startsWith(ctx.m_pepSeq,cellTemp{2})
				continue;
			else
				fixedPosMod{1,1} = 1;
				fixedPosMod{1,2} = ctx.m_fixedModNameMass{idx,1};
				fixedPosMod{1,3} = ctx.m_fixedModNameMass{idx,3};
			end
		else
			fixedPosMod{1,1} = 0;
			fixedPosMod{1,2} = ctx.m_fixedModNameMass{idx,1};
			fixedPosMod{1,3} = ctx.m_fixedModNameMass{idx,3};
		end
	elseif contains(strSpecfct,'C-term')
		cellTemp=split(strSpecfct,'C-term');
		if isequal(cellTemp{1},'Protein') && ~ctx.m_isProtC
			continue;
		end
		if ~isempty(cellTemp{2})
			if ~endsWith(ctx.m_pepSeq,cellTemp{2})
				continue;   % C-terminal specificity, with amino acid restriction
			else
				fixedPosMod{length(ctx.m_pepSeq)+2,1} = length(ctx.m_pepSeq);
				fixedPosMod{length(ctx.m_pepSeq)+2,2} = ctx.m_fixedModNameMass{idx,1};
				fixedPosMod{length(ctx.m_pepSeq)+2,3} = ctx.m_fixedModNameMass{idx,3};
			end
		else
			fixedPosMod{length(ctx.m_pepSeq)+2,1} = length(ctx.m_pepSeq)+1;
			fixedPosMod{length(ctx.m_pepSeq)+2,2} = ctx.m_fixedModNameMass{idx,1};
			fixedPosMod{length(ctx.m_pepSeq)+2,3} = ctx.m_fixedModNameMass{idx,3};
		end
	else
		% Only one amino acid
		vecPos=strfind(ctx.m_pepSeq,strSpecfct);
		for iPos=1:length(vecPos)
			fixedPosMod{vecPos(iPos)+1,1} = vecPos(iPos);
			fixedPosMod{vecPos(iPos)+1,2} = ctx.m_fixedModNameMass{idx,1};
			fixedPosMod{vecPos(iPos)+1,3} = ctx.m_fixedModNameMass{idx,3};
		end
	end
end
fixedPosMod(cellfun(@isempty,fixedPosMod(:,1)),:)=[];
end