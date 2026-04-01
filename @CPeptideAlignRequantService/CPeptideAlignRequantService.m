classdef CPeptideAlignRequantService < handle
    % Service for peptide-level alignment and re-quantification.

    properties (Access = private)
        m_cfg
        m_fixedModNameMass
        m_variableModNameMass
        m_cMgfDatasetIO
        m_cMs12DatasetIO
        m_cMsFileMapper
        m_pepProtService
    end

    methods
        function obj = CPeptideAlignRequantService(cfg)
            % Input:
            %   cfg (struct)
            %       minimal config struct
            if ~isstruct(cfg)
                error('CPeptideAlignRequantService:InvalidConstructorArgs', ...
                    'Expected a config struct.');
            end
            obj.m_cfg = cfg;

            [obj.m_fixedModNameMass, obj.m_variableModNameMass] = CModificationRegistry.fromConfig(obj.m_cfg);
            obj.m_cMgfDatasetIO = CMgfDatasetIO(obj.m_cfg.spec_dir_path);
            obj.m_cMs12DatasetIO = CMS12DatasetIO(obj.m_cfg.spec_dir_path, obj.m_cfg.ms1_tolerance);
            obj.m_cMsFileMapper = CMsFileMapper(obj.m_cfg.spec_dir_path);
            obj.m_pepProtService = CPepProtService(obj.m_cfg.fasta_file_path, ...
                obj.m_cfg.regular_express, obj.m_cfg.filtered_res_file_path);
        end

        function run(obj)
            % Run peptide-level alignment and re-quantification.
            % Alignment strategy and options are provided via cfg fields:
            %   - cfg.align_strategy_obj
            %   - cfg.align_options
            % Recognized align_options keys:
            %       - peptide_quant_res_path (1 x 1 char/string)
            %           initial peptide quant result path (deprecated; use cfg.peptide_quant_res_path)
            %       - min_psm (1 x 1 double)
            %           minimum anchors/PSMs for alignment (passed to aligner)
            %       - num_bins (1 x 1 double)
            %           number of RT bins in local transform fitting (passed to aligner)
            %       - min_per_bin (1 x 1 double)
            %           minimum anchors per RT bin (passed to aligner)
            %       - outlier_k (1 x 1 double)
            %           outlier threshold factor for transform fitting (passed to aligner)
            %       - outlier_method (1 x 1 char/string)
            %           outlier method, e.g. 'mad'/'std' (passed to aligner)
            %       - rt_sigma (1 x 1 double)
            %           RT Gaussian sigma for peak selection (passed to aligner)
            %       - max_rt_residual (1 x 1 double)
            %           max RT residual for candidate pairing (passed to aligner)
            %       - dead_time_min (1 x 1 double)
            %           minimum allowed RT start (passed to aligner)
            %       - unaligned_imp_action (1 x 1 char/string)
            %           behavior when no valid aligned pair is selected:
            %             'drop' -> remove IMP from pep_rtrange_map (default)
            %             'zero-label' -> keep peaks but set all check_label = 0
            %       unknown keys are ignored by this service and forwarded downstream
            cfg = obj.m_cfg;
            align_strategy = cfg.align_strategy_obj;
            align_options = cfg.align_options;

            if isempty(align_strategy)
                error('CPeptideAlignRequantService:MissingAlignStrategy', ...
                    'align_strategy is required.');
            end
            if isempty(align_options)
                error('CPeptideAlignRequantService:MissingAlignOptions', ...
                    'align_options is required.');
            end

            if isempty(cfg.filtered_res_file_path)
                error('CPeptideAlignRequantService:MissingFilteredResult', ...
                    'FDR filtered result path is required for alignment.');
            end

            msms_res_path = cfg.msms_res_path;
            msms_result = CMS2ResultIO.read(msms_res_path);

            quant_output_path = cfg.peptide_quant_res_path;
            obj.ensureParentDirForFile(quant_output_path);

            proc_cfg = CIMPProcessingExecutorConfig(struct( ...
                'ms12DatasetIO', obj.m_cMs12DatasetIO, ...
                'ms1_tolerance', cfg.ms1_tolerance, ...
                'minMSMSnum', cfg.min_MSMS_num, ...
                'alpha', cfg.alpha, ...
                'resFilterThres', 0));  % cfg.result_filter_threshold
            proc_executor = CIMPProcessingExecutor(proc_cfg);

            rawIdentManagers = cell(1, length(msms_result.Peptides));
            base_groups_by_peptide = cell(1, length(msms_result.Peptides));
            prot_name_pos_by_peptide = cell(1, length(msms_result.Peptides));
            print_progress = CPrintProgress(length(msms_result.Peptides), 'build_raw_ident_for_alignment');
            for i_pep = 1:length(msms_result.Peptides)
                print_progress = print_progress.update_show(i_pep);
                rawIdentManager = obj.buildRawIdentManagerFromSpectrumList(msms_result.Peptides(i_pep).spectrum_list);
                rawIdentManagers{i_pep} = rawIdentManager;
                base_groups_by_peptide{i_pep} = proc_executor.buildBaseGroups(rawIdentManager);

                peptide_sequence = msms_result.Peptides(i_pep).peptide_sequence;
                prot_name_pos_by_peptide{i_pep} = obj.m_pepProtService.get_protein_name_pos(peptide_sequence);
            end
            print_progress.last_update();

            quant_report = CIMPQuantReport();
            print_progress = CPrintProgress(length(msms_result.Peptides), 'peptide_quant_before_alignment');
            CLogger.info('Quantifying at peptide level before alignment...');
            for i_pep = 1:length(msms_result.Peptides)
                print_progress = print_progress.update_show(i_pep);
                cell_prot_name_pos = prot_name_pos_by_peptide{i_pep};
                rawIdentManager = rawIdentManagers{i_pep};
                block = proc_executor.quantifyPeptideBlock(...
                    cell_prot_name_pos, rawIdentManager, base_groups_by_peptide{i_pep});
                quant_report = quant_report.append_block(block);
            end
            print_progress.last_update();
            CLogger.info('Peptide-level quantification before alignment done.');

            CIMPQuantResultIO.write(quant_report, quant_output_path);
            % quant_report = CIMPQuantResultIO.read(quant_output_path);
            base_pep_rtrange_map = quant_report.build_pep_rtrange_map();

            anchor_selector = CXICAlignAnchorSelector();
            aligner = CXICAligner(anchor_selector, align_options);
            align_cfg = CIMPXICAlignRequantExecutorConfig(struct( ...
                'ms12DatasetIO', obj.m_cMs12DatasetIO, ...
                'ms1_tolerance', cfg.ms1_tolerance, ...
                'minMSMSnum', cfg.min_MSMS_num, ...
                'resFilterThres', cfg.result_filter_threshold, ...
                'aligner', aligner, ...
                'align_strategy', align_strategy, ...
                'align_options', align_options));
            align_executor = CIMPXICAlignRequantExecutor(align_cfg);

            [pep_rtrange_map, align_report] = align_executor.buildAlignedRtRangeMap( ...
                cfg.filtered_res_file_path, rawIdentManagers, base_pep_rtrange_map, base_groups_by_peptide);

            obj.ensureParentDirForFile(cfg.alignment_report_path);
            align_executor.writeAlignmentReport(align_report, cfg.alignment_report_path);
            
            if ~isempty(cfg.align_requant_rt_stats_path)
                obj.ensureParentDirForFile(cfg.align_requant_rt_stats_path);
                CIMPQuantStats.rt_sorted_stats('init');
                stats_cleanup = onCleanup(@() CIMPQuantStats.rt_sorted_stats('flush', cfg.align_requant_rt_stats_path)); 
            end

            output_path = cfg.requant_output_path;
            obj.ensureParentDirForFile(output_path);
            report = CIMPQuantReport();
            print_progress = CPrintProgress(length(msms_result.Peptides), 'peptide_requant_aligned');

            CLogger.info('Re-quantifying at peptide level (aligned)...');
            for i_pep = 1:length(msms_result.Peptides)
                print_progress = print_progress.update_show(i_pep);
                cell_prot_name_pos = prot_name_pos_by_peptide{i_pep};
                rawIdentManager = rawIdentManagers{i_pep};
                block = proc_executor.requantifyPeptideBlock(...
                    cell_prot_name_pos, rawIdentManager, pep_rtrange_map, base_groups_by_peptide{i_pep});
                report = report.append_block(block);
            end
            print_progress.last_update();
            CLogger.info('Aligned peptide-level re-quantification done.');

            CIMPQuantResultIO.write(report, output_path);
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


        function ensureParentDirForFile(~, file_path)
            % Ensure parent directory exists for a target file path.
            % Input:
            %   file_path (1 x 1 char/string)
            %       target file path
            % Output:
            %   (none)
            if isempty(file_path)
                return;
            end

            parent_dir = fileparts(char(string(file_path)));
            if isempty(parent_dir)
                return;
            end
            CPathResolver.ensureDir(parent_dir);
        end

    end
end
