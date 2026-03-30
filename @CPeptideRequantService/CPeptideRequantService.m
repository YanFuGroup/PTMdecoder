classdef CPeptideRequantService < handle
    % Stage service for peptide-level re-quantification.

    properties (Access = private)
        m_msms_cfg
    end

    methods
        function obj = CPeptideRequantService(msms_cfg)
            % Input:
            %   msms_cfg (struct)
            %       config for peptide-level re-quant processor
            if nargin < 1 || isempty(msms_cfg)
                error('CPeptideRequantService:MissingConfig', ...
                    'msms_cfg must be provided.');
            end
            obj.m_msms_cfg = msms_cfg;
        end

        function run(obj)
            % Run peptide-level re-quantification stage.
            cfg = obj.m_msms_cfg;

            [fixedModNameMass, variableModNameMass] = CModificationRegistry.fromConfig(cfg);

            checked_pep_path = CPathResolver.resolveFilePath(cfg.output_dir_path, ...
                'report_peptide_all_checked.txt', cfg.checked_peptides_res_path);

            report = CIMPQuantResultIO.read(checked_pep_path);
            pep_rtrange_map = report.build_pep_rtrange_map();
            if isempty(pep_rtrange_map)
                error(['The checked peptide result file "', checked_pep_path, '" is empty!']);
            end

            CPathResolver.ensureDir(cfg.output_dir_path);

            each_PSM_results_path = CPathResolver.resolveFilePath(cfg.output_dir_path, ...
                'report_msms.txt', cfg.msms_res_path);
            output_path = CPathResolver.resolveFilePath(cfg.output_dir_path, ...
                'report_peptide_all_requant.txt', '');

            pepProtService = CPepProtService(cfg.fasta_file_path, cfg.regular_express, cfg.filtered_res_file_path);
            cMgfDatasetIO = CMgfDatasetIO(cfg.spec_dir_path);
            cMs12DatasetIO = CMS12DatasetIO(cfg.spec_dir_path, cfg.ms1_tolerance);
            cMsFileMapper = CMsFileMapper(cfg.spec_dir_path);

            msms_result = CMS2ResultIO.read(each_PSM_results_path);
            print_progress = CPrintProgress(length(msms_result.Peptides), 'peptide_requantification');
            pipeline_cfg = CIMPProcessingExecutorConfig(struct( ...
                'ms12DatasetIO', cMs12DatasetIO, ...
                'ms1_tolerance', cfg.ms1_tolerance, ...
                'minMSMSnum', cfg.min_MSMS_num, ...
                'alpha', cfg.alpha, ...
                'resFilterThres', cfg.result_filter_threshold));
            executor = CIMPProcessingExecutor(pipeline_cfg);
            CIMPQuantStats.rt_sorted_stats('init');
            stats_cleanup = onCleanup(@() CIMPQuantStats.rt_sorted_stats('flush', ...
                CPathResolver.resolveFilePath(cfg.output_dir_path, 'requant_rt_sorted_stats.mat', '')));

            deps = CPeptideRawIdentAssembler.createSpectrumListDeps( ...
                cMgfDatasetIO, cMs12DatasetIO, cMsFileMapper, cfg.ms1_tolerance, ...
                fixedModNameMass, variableModNameMass, cfg.msms_stability_filter);

            report_requant = CIMPQuantReport();
            CLogger.info('Re-quantifying at peptide level...');
            for idx_psf = 1:length(msms_result.Peptides)
                print_progress = print_progress.update_show(idx_psf);
                peptide_sequence = msms_result.Peptides(idx_psf).peptide_sequence;
                cell_prot_name_pos = pepProtService.get_protein_name_pos(peptide_sequence);
                rawIdentManager = CPeptideRawIdentAssembler.buildFromSpectrumList( ...
                    msms_result.Peptides(idx_psf).spectrum_list, deps);
                block = executor.requantifyPeptideBlock(cell_prot_name_pos, rawIdentManager, pep_rtrange_map);
                report_requant = report_requant.append_block(block);
            end
            print_progress.last_update();
            CLogger.info('Peptide-level re-quantification done.');

            CIMPQuantResultIO.write(report_requant, output_path);
        end
    end
end
