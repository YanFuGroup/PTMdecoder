classdef CModificationRegistry
    % Centralized helper for modification map loading and mod-string parsing.

    methods (Static)
        function [fixedModNameMass, variableModNameMass, mapModification] = fromConfig(cfg)
            % Build fixed/variable modification tables from stage config.
            % Input:
            %   cfg (struct)
            %       must contain: mod_file_path, fixed_mod, variable_mod
            % Output:
            %   fixedModNameMass (Nf x 3 cell)
            %   variableModNameMass (Nv x 3 cell)
            %   mapModification (containers.Map)
            required_fields = {'mod_file_path', 'fixed_mod', 'variable_mod'};
            for i = 1:numel(required_fields)
                if ~isfield(cfg, required_fields{i})
                    error('CModificationRegistry:MissingConfigField', ...
                        'Missing required config field: %s', required_fields{i});
                end
            end

            mapModification = CModificationRegistry.loadMap(cfg.mod_file_path);
            fixedModNameMass = CModificationRegistry.buildModNameMass(cfg.fixed_mod, mapModification);
            variableModNameMass = CModificationRegistry.buildModNameMass(cfg.variable_mod, mapModification);
        end

        function mapModification = loadMap(strPathname)
            % Read pFind style modification ini and build name->mono-mass map.
            mapModification = containers.Map();

            fid = fopen(strPathname,'r');
            if 0>=fid
                error('CModificationRegistry:OpenFileFailed', ...
                    'Failed to open modify.ini file!');
            end
            file_cleanup = onCleanup(@() fclose(fid));

            % Skip the first line of header, and read the number of modifications from the second line.
            fgetl(fid);
            strLine = fgetl(fid);
            [~,strEnd] = strtok(strLine,'=');
            nModifications = str2double(strEnd(2:end));

            nModify = 0;
            while ~feof(fid)
                strLine = fgetl(fid);
                if isempty(strLine)
                    continue;
                end

                if strncmpi(strLine,'name',4)
                    % If it is a 'name=' line, then read the next line, which is its detailed information line
                    strLine = fgetl(fid);
                    [strName,strValue] = strtok(strLine,'=');

                    % The mono mass of the modification is between the second and third blank space ' ', convert it to a numeric type
                    Idx = strfind(strValue,' ');
                    strValue = str2double(strValue(Idx(2)+1:Idx(3)-1));

                    mapModification(strName) = strValue;
                    nModify = nModify + 1;
                end
            end

            if nModifications ~= nModify
                error('CModificationRegistry:InvalidHeaderCount', ...
                    'Failed to read modify.ini file! The number of read modification is not consistent with the header.');
            end

            mapModification('null') = 0;
            mapModification('NULL') = 0;
            mapModification('Null') = 0;
            mapModification('unknown') = 0;

            clear file_cleanup;
        end

        function modNameMass = buildModNameMass(modificationTypes, mapModification)
            % Build modification name/mass rows from setting string.
            % Input:
            %   modificationTypes (1 x 1 char/string)
            %       semicolon-separated modification declarations
            %   mapModification (containers.Map)
            %       modification declaration to mass map
            % Output:
            %   modNameMass (N x 3 cell)
            %       {mod_name, specificity, mass}
            modNameMass = [];
            if isempty(modificationTypes)
                return
            end

            S_modificationTypes = regexp(modificationTypes,';','split');
            modNameMass = cell(length(S_modificationTypes),3);
            for i = 1:length(S_modificationTypes)
                if isempty(S_modificationTypes{i})
                    continue;
                end

                left_brac_pos = strfind(S_modificationTypes{i},'[');
                if length(left_brac_pos)>2
                    error(['Unexpected modification: ',S_modificationTypes{i}, ...
                        'The modification string are expected to be in either ' ...
                        '"Carbamidomethyl[C]" or "ICPL_13C(6)[K](NIC_13C(6)[K])"']);
                elseif length(left_brac_pos)==2
                    right_brac_pos = strfind(S_modificationTypes{i},']');
                    modNameMass{i,1} = [S_modificationTypes{i}(1:left_brac_pos(1)-1),...
                        S_modificationTypes{i}(right_brac_pos(1)+1:left_brac_pos(2)),')'];
                    modNameMass{i,2} = S_modificationTypes{i}(left_brac_pos(1)+1:...
                        right_brac_pos(1)-1);
                else
                    modNameMass{i,1} = S_modificationTypes{i}(1:left_brac_pos(1)-1);
                    modNameMass{i,2} = S_modificationTypes{i}(left_brac_pos(1)+1:end-1);
                end

                modNameMass{i,3} = mapModification(S_modificationTypes{i});
            end
            modNameMass(cellfun(@isempty,modNameMass(:,1)),:) = [];
        end
    end
end
