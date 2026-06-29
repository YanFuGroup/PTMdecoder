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
            [obj.m_fixedModNameMass, obj.m_variableModNameMass] = CModificationRegistry.fromConfig(obj.m_cfg);
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
            if isfield(cfg, 'xic_layout')
                xic_layout = cfg.xic_layout;
            else
                xic_layout = [];
            end

            checked_pep_path = CPathResolver.resolveFilePath(cfg.output_dir_path, ...
                'report_peptide_all_checked.txt', cfg.checked_peptides_res_path);
            report = CIMPQuantResultIO.read(checked_pep_path);
            pep_rtrange_map = report.build_pep_rtrange_map();
            if isempty(pep_rtrange_map)
                CLogger.warn('The checked peptide result file "%s" is empty!', checked_pep_path);
            end

            CPathResolver.ensureDir(dir_save);

            each_PSM_results_path = CPathResolver.resolveFilePath(cfg.output_dir_path, ...
                'report_msms.txt', cfg.msms_res_path);
            msms_result = CMS2ResultIO.read(each_PSM_results_path);

            print_progress = CPrintProgress(length(msms_result.Peptides), 'draw_xic');
            pipeline_cfg = CIMPProcessingExecutorConfig(struct( ...
                'ms12DatasetIO', obj.m_cMs12DatasetIO, ...
                'ms1_tolerance', cfg.ms1_tolerance, ...
                'minMSMSnum', cfg.min_MSMS_num, ...
                'alpha', cfg.alpha, ...
                'resFilterThres', cfg.result_filter_threshold));
            executor = CIMPProcessingExecutor(pipeline_cfg);

            CLogger.info('Drawing XIC...');
            for idx_psf = 1:length(msms_result.Peptides)
                print_progress = print_progress.update_show(idx_psf);
                rawIdentManager = obj.buildRawIdentManagerFromSpectrumList(msms_result.Peptides(idx_psf).spectrum_list);
                executor.drawImpXicForBlock(rawIdentManager, pep_rtrange_map, ...
                    dir_save, color_map, legend_map, xic_layout);
            end
            print_progress.last_update();
            CLogger.info('Drawing XIC done.');
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
            deps = CPeptideRawIdentAssembler.createSpectrumListDeps( ...
                obj.m_cMgfDatasetIO, obj.m_cMs12DatasetIO, obj.m_cMsFileMapper, obj.m_cfg.ms1_tolerance, ...
                obj.m_fixedModNameMass, obj.m_variableModNameMass, obj.m_cfg.msms_stability_filter);
            rawIdentManager = CPeptideRawIdentAssembler.buildFromSpectrumList(spectrum_list, deps);
        end

    end
end
