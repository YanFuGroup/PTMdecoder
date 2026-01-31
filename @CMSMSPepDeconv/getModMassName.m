function [modNameMass]=getModMassName(~,modificationTypes,mapModification)
% Use the modification types input by the user to search the modification library and create a matrix of modification names corresponding to modification masses
% Input:
%   modificationTypes (1 x 1 char/string)
%       user-specified modification types, separated by semicolons
%       e.g., 'Carbamidomethyl[C];Oxidation[M]'
%   mapModification (containers.Map)
%       modification library map
% Output:
%   modNameMass (M x 3 cell)
%       modification name, specificity, and mass
modNameMass=[];
if isempty(modificationTypes)
    return
end
S_modificationTypes=regexp(modificationTypes,';','split');
% The first element of each row is the modification name, the second is the specific site, and the third is the modification mass
modNameMass=cell(length(S_modificationTypes),3);
for i=1:length(S_modificationTypes)
    if isempty(S_modificationTypes{i})
        continue;
    end
    left_brac_pos = strfind(S_modificationTypes{i},'[');
    if length(left_brac_pos)>2
        error(['Unexpected modification: ',S_modificationTypes{i}, ...
            'The modification string are expected to be in either ' ...
            '"Carbamidomethyl[C]" or "ICPL_13C(6)[K](NIC_13C(6)[K])"']);
    elseif length(left_brac_pos)==2
        %　some modification may be like "ICPL_13C(6)[K](NIC_13C(6)[K])"
        right_brac_pos = strfind(S_modificationTypes{i},']');
        modNameMass{i,1} = [S_modificationTypes{i}(1:left_brac_pos(1)-1),...
            S_modificationTypes{i}(right_brac_pos(1)+1:left_brac_pos(2)),')'];
        modNameMass{i,2} = S_modificationTypes{i}(left_brac_pos(1)+1:...
            right_brac_pos(1)-1);
    else
        % normal format like "Acetyl[S]"
        modNameMass{i,1} = S_modificationTypes{i}(1:left_brac_pos(1)-1);
        modNameMass{i,2} = S_modificationTypes{i}(left_brac_pos(1)+1:end-1);
    end
	modNameMass{i,3} = mapModification(S_modificationTypes{i});
end
modNameMass(cellfun(@isempty,modNameMass(:,1)),:)=[];
end
