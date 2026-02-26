classdef CXICDrawService < handle
    % Service for drawing XIC from checked peptide RT ranges.

    properties (Access = private)
        m_cfg
        m_fixedModNameMass
        m_variableModNameMass
        m_cMgfDatasetIO
        m_cMs12DatasetIO
        m_cMsFileMapper
    end

    methods
        function obj = CXICDrawService(cfg)
            % Input:
            %   cfg (struct)
            %       minimal config struct
            if ~isstruct(cfg)
                error('CXICDrawService:InvalidConstructorArgs', ...
                    'Expected a config struct.');
            end
            obj.m_cfg = cfg;
            mapModification = readModifyInfo(obj.m_cfg.mod_file_path);
            obj.m_fixedModNameMass = obj.getModMassName(obj.m_cfg.fixed_mod, mapModification);
            obj.m_variableModNameMass = obj.getModMassName(obj.m_cfg.variable_mod, mapModification);
            obj.m_cMgfDatasetIO = CMgfDatasetIO(obj.m_cfg.spec_dir_path);
            obj.m_cMs12DatasetIO = CMS12DatasetIO(obj.m_cfg.spec_dir_path, obj.m_cfg.ms1_tolerance);
            obj.m_cMsFileMapper = CMsFileMapper(obj.m_cfg.spec_dir_path);
        end

        function run(obj, dir_save, color_map, legend_map)
            % Draw XIC figures from checked peptide RT ranges.
            % Input:
            %   dir_save (1 x 1 char/string)
            %       output directory for figures
            %   color_map (containers.Map, optional)
            %       map from modified peptide string to RGB color (1 x 3 double)
            %   legend_map (containers.Map, optional)
            %       map from modified peptide string to legend display text
            %       if empty, default plotting color/legend behavior is used
            if nargin < 4
                legend_map = [];
            end
            if nargin < 3
                color_map = [];
            end

            cfg = obj.m_cfg;

            checked_pep_path = CPathResolver.resolveFilePath(cfg.output_dir_path, ...
                'report_peptide_all_checked.txt', cfg.checked_peptides_res_path);
            report = CIMPQuantResultIO.read(checked_pep_path);
            pep_rtrange_map = report.build_pep_rtrange_map();
            if isempty(pep_rtrange_map)
                warning(['The checked peptide result file "', checked_pep_path, '" is empty!']);
            end

            CPathResolver.ensureDir(dir_save);

            each_PSM_results_path = CPathResolver.resolveFilePath(cfg.output_dir_path, ...
                'report_msms.txt', cfg.msms_res_path);
            msms_result = CMS2ResultIO.read(each_PSM_results_path);

            print_progress = CPrintProgress(length(msms_result.Peptides));
            pipeline_cfg = CIMPProcessingExecutorConfig(struct( ...
                'ms12DatasetIO', obj.m_cMs12DatasetIO, ...
                'ms1_tolerance', cfg.ms1_tolerance, ...
                'minMSMSnum', cfg.min_MSMS_num, ...
                'alpha', cfg.alpha, ...
                'resFilterThres', cfg.result_filter_threshold));
            executor = CIMPProcessingExecutor(pipeline_cfg);
            pipeline = CPeptideLevelPipeline(executor);

            fprintf('Drawing XIC...');
            for idx_psf = 1:length(msms_result.Peptides)
                print_progress = print_progress.update_show(idx_psf);
                rawIdentManager = obj.buildRawIdentManagerFromSpectrumList(msms_result.Peptides(idx_psf).spectrum_list);
                pipeline.drawImpXicForBlock(rawIdentManager, pep_rtrange_map, dir_save, color_map, legend_map);
            end
            print_progress.last_update();
            fprintf('done.\n');
        end
    end

    methods (Access = private)
        function rawIdentManager = buildRawIdentManagerFromSpectrumList(obj, spectrum_list)
            % Build raw identification manager from spectrum list.
            % Input:
            %   spectrum_list (struct array)
            %       spectrum list from MSMS result
            % Output:
            %   rawIdentManager (CIMPRawIdentManager)
            %       assembled raw identification manager
            deps = struct( ...
                'getProfilesFunc', @(dataset_name, spectrum_name) obj.getProfiles(dataset_name, spectrum_name), ...
                'fixedModNameMass', {obj.m_fixedModNameMass}, ...
                'variableModNameMass', {obj.m_variableModNameMass});
            rawIdentManager = CPeptideRawIdentAssembler.buildFromSpectrumList(spectrum_list, deps);
        end

        function [isorts,c_ref_isointens,cur_mz,cur_ch] = getProfiles(obj, mgf_name, spectrum_name)
            % Read MS1 profile around one MS2 spectrum precursor.
            % Input:
            %   mgf_name (1 x 1 char/string)
            %       dataset file name in MGF
            %   spectrum_name (1 x 1 char/string)
            %       spectrum name in MGF
            % Output:
            %   isorts (1 x 1 double)
            %       MS1 retention time of matched precursor scan
            %   c_ref_isointens (1 x 1 double)
            %       reference isotope intensity near precursor m/z
            %   cur_mz (1 x 1 double)
            %       precursor m/z from MGF
            %   cur_ch (1 x 1 double/int)
            %       precursor charge from MGF
            spec_name = regexp(spectrum_name,'\.','split');
            MS2ScanI = str2double(spec_name{2});
            [~, cur_ch, cur_mz] = obj.m_cMgfDatasetIO.read_oneSpec(mgf_name,spectrum_name);

            mgf_stem = erase(mgf_name,'.mgf');
            ms12_stem = obj.m_cMsFileMapper.get_ms1_stem(mgf_stem);

            MS2_index = obj.m_cMs12DatasetIO.m_mapNameMS2Index(ms12_stem);
            idx_cur_scan = MS2_index(:,2)==MS2ScanI;
            MS1Scan = MS2_index(idx_cur_scan,1);
            MS1_index = obj.m_cMs12DatasetIO.m_mapNameMS1Index(ms12_stem);
            MS1_peaks = obj.m_cMs12DatasetIO.m_mapNameMS1Peaks(ms12_stem);
            index_starts_MS1 = [1;MS1_index(1:size(MS1_index,1),3)];
            ino = find(MS1_index(:,1)==MS1Scan);
            isorts = MS1_index(ino,2);
            IX = index_starts_MS1(ino):index_starts_MS1(ino+1)-1;
            mz = MS1_peaks(IX,1);
            inten = MS1_peaks(IX,2);
            if obj.m_cfg.ms1_tolerance.isppm
                ptol = obj.m_cfg.ms1_tolerance.value*cur_mz*1e-6;
            else
                ptol = obj.m_cfg.ms1_tolerance.value;
            end
            c_ptol = min([ptol,0.3]);
            c_ref_isointens = max(inten(abs(mz-cur_mz)<c_ptol));
            if isempty(c_ref_isointens)
                c_ref_isointens = 0;
            end
        end

        function modNameMass = getModMassName(~, modificationTypes, mapModification)
            % Build modification name/mass mapping list from setting string.
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
