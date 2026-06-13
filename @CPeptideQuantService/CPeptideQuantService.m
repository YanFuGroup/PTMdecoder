classdef CPeptideQuantService < handle
    % Stage service for peptide-level quantification.

    properties (Access = private)
        m_msms_cfg
    end

    methods
        function obj = CPeptideQuantService(msms_cfg)
            % Input:
            %   msms_cfg (struct)
            %       config for peptide-level processor
            if nargin < 1 || isempty(msms_cfg)
                error('CPeptideQuantService:MissingConfig', ...
                    'msms_cfg must be provided.');
            end
            obj.m_msms_cfg = msms_cfg;
        end

        function run(obj)
            % Run peptide-level quantification stage.
            cfg = obj.m_msms_cfg;

            [fixedModNameMass, variableModNameMass] = CModificationRegistry.fromConfig(cfg);

            pepProtService = CPepProtService(cfg.fasta_file_path, cfg.regular_express, cfg.filtered_res_file_path);
            cMgfDatasetIO = CMgfDatasetIO(cfg.spec_dir_path);
            cMs12DatasetIO = CMS12DatasetIO(cfg.spec_dir_path, cfg.ms1_tolerance);
            cMsFileMapper = CMsFileMapper(cfg.spec_dir_path);

            CPathResolver.ensureDir(cfg.output_dir_path);

            each_PSM_results_path = CPathResolver.resolveFilePath(cfg.output_dir_path, ...
                'report_msms.txt', cfg.msms_res_path);
            % TODO: Set the output path for peptide level results
            each_peptide_results_path = CPathResolver.resolveFilePath(cfg.output_dir_path, ...
                'report_peptide_all.txt', '');

            report = CIMPQuantReport();
            msms_result = CMS2ResultIO.read(each_PSM_results_path);
            print_progress = CPrintProgress(length(msms_result.Peptides), 'IMP quantification');

            pipeline_cfg = CIMPProcessingExecutorConfig(struct( ...
                'ms12DatasetIO', cMs12DatasetIO, ...
                'ms1_tolerance', cfg.ms1_tolerance, ...
                'minMSMSnum', cfg.min_MSMS_num, ...
                'alpha', cfg.alpha, ...
                'resFilterThres', cfg.result_filter_threshold, ...
                'minXicNonzeroPoints', cfg.min_xic_nonzero_points));
            executor = CIMPProcessingExecutor(pipeline_cfg);
            CIMPQuantStats.rt_sorted_stats('init');
            CIMPQuantStats.quant_group_stats('init', []);
            stats_cleanup = onCleanup(@() CIMPQuantStats.rt_sorted_stats('flush', ...
                CPathResolver.resolveFilePath(cfg.output_dir_path, 'rt_sorted_stats.mat', '')));

            deps = CPeptideRawIdentAssembler.createSpectrumListDeps( ...
                cMgfDatasetIO, cMs12DatasetIO, cMsFileMapper, cfg.ms1_tolerance, ...
                fixedModNameMass, variableModNameMass, cfg.msms_stability_filter);

            CLogger.info('Quantifying at peptide level...');
            total_filter_stats = struct( ...
                'total_spectra', 0, ...
                'kept_spectra', 0, ...
                'dropped_spectra_jaccard', 0, ...
                'total_imp_candidates', 0, ...
                'kept_imp_candidates', 0, ...
                'dropped_imp_candidates', 0);
            for idx_psf = 1:length(msms_result.Peptides)
                print_progress = print_progress.update_show(idx_psf);
                peptide_sequence = msms_result.Peptides(idx_psf).peptide_sequence;
                cell_prot_name_pos = pepProtService.get_protein_name_pos(peptide_sequence);
                [rawIdentManager, filter_stats] = CPeptideRawIdentAssembler.buildFromSpectrumList( ...
                    msms_result.Peptides(idx_psf).spectrum_list, deps);
                total_filter_stats.total_spectra = total_filter_stats.total_spectra + filter_stats.total_spectra;
                total_filter_stats.kept_spectra = total_filter_stats.kept_spectra + filter_stats.kept_spectra;
                total_filter_stats.dropped_spectra_jaccard = total_filter_stats.dropped_spectra_jaccard + filter_stats.dropped_spectra_jaccard;
                total_filter_stats.total_imp_candidates = total_filter_stats.total_imp_candidates + filter_stats.total_imp_candidates;
                total_filter_stats.kept_imp_candidates = total_filter_stats.kept_imp_candidates + filter_stats.kept_imp_candidates;
                total_filter_stats.dropped_imp_candidates = total_filter_stats.dropped_imp_candidates + filter_stats.dropped_imp_candidates;
                block = executor.quantifyPeptideBlock(cell_prot_name_pos, rawIdentManager);
                report = report.append_block(block);
            end
            print_progress.last_update();
            CIMPQuantStats.log_quant_group_summary('Peptide-level');

            if cfg.msms_stability_filter.enabled
                CLogger.info(['MSMS stability filter summary: spectra(total=%d, kept=%d, dropped_by_jaccard=%d), ', ...
                    'imps(total=%d, kept=%d, dropped=%d).'], ...
                    total_filter_stats.total_spectra, total_filter_stats.kept_spectra, total_filter_stats.dropped_spectra_jaccard, ...
                    total_filter_stats.total_imp_candidates, total_filter_stats.kept_imp_candidates, total_filter_stats.dropped_imp_candidates);
            end

            CLogger.info('Peptide-level quantification done.');

            CIMPQuantResultIO.write(report, each_peptide_results_path);
        end
    end
end
