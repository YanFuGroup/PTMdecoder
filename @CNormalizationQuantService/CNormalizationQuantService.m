classdef CNormalizationQuantService < handle
    % Stage service for normalization-peptide quantification.

    properties (Access = private)
        m_cfg
    end

    methods
        function obj = CNormalizationQuantService(cfg)
            % Input:
            %   cfg (struct)
            %       - msms_cfg (struct, minimal normalization config)
            %           required fields:
            %           spec_dir_path, ms1_tolerance, alpha,
            %           result_filter_threshold, output_dir_path
            %       - peptide_list (K x 1 cell)
            %       - prot_list (K x 1 cell)
            %       - filtered_res_file_path (1 x 1 char/string)
            %       - output_file_name (1 x 1 char/string)
            if nargin < 1 || isempty(cfg)
                error('CNormalizationQuantService:MissingConfig', ...
                    'cfg must be provided.');
            end
            obj.m_cfg = cfg;
        end

        function run(obj)
            % Run normalization peptide quantification stage.
            cfg = obj.m_cfg;
            msms_cfg = cfg.msms_cfg;

            peptide_list = cfg.peptide_list;
            prot_list = cfg.prot_list;
            if length(peptide_list) ~= length(prot_list)
                error('CNormalizationQuantService:InvalidInput', ...
                    'peptide_list and prot_list must have the same length');
            end

            CPathResolver.ensureDir(msms_cfg.output_dir_path);

            cMs12DatasetIO = CMS12DatasetIO(msms_cfg.spec_dir_path, msms_cfg.ms1_tolerance);

            pep_quant = cell(length(peptide_list), 1);
            for i_list = 1:length(peptide_list)
                pep_quant{i_list} = CIMPRawIdentManager();
            end

            FDRfilteredResults = CFdrFilteredResultIO.read(cfg.filtered_res_file_path);
            entries = FDRfilteredResults.entries;
            deps = struct( ...
                'ms12DatasetIO', cMs12DatasetIO, ...
                'ms1_tolerance', msms_cfg.ms1_tolerance);

            fprintf('Reading %s...', msms_cfg.spec_dir_path);
            pep_quant = CPeptideRawIdentAssembler.buildFromFdrEntries(entries, peptide_list, pep_quant, deps);
            fprintf('done.\n');

            pipeline_cfg = CIMPProcessingExecutorConfig(struct( ...
                'ms12DatasetIO', cMs12DatasetIO, ...
                'ms1_tolerance', msms_cfg.ms1_tolerance, ...
                'minMSMSnum', 1, ...
                'alpha', msms_cfg.alpha, ...
                'resFilterThres', msms_cfg.result_filter_threshold));
            executor = CIMPProcessingExecutor(pipeline_cfg);
            CIMPQuantStats.rt_sorted_stats('init');
            stats_cleanup = onCleanup(@() CIMPQuantStats.rt_sorted_stats('flush', ...
                CPathResolver.resolveFilePath(msms_cfg.output_dir_path, 'rt_sorted_stats.mat', '')));
            report = CIMPQuantReport();

            fprintf('Quantifying %s...', msms_cfg.spec_dir_path);
            for i_list = 1:length(peptide_list)
                block = executor.quantifyPeptideBlock({prot_list{i_list}, -1}, pep_quant{i_list});
                report = report.append_block(block);
            end
            fprintf('done.\n');

            % TODO: Set the output path for peptide level results
            output_path = CPathResolver.resolveFilePath(msms_cfg.output_dir_path, cfg.output_file_name, '');
            CIMPQuantResultIO.write(report, output_path);
        end
    end
end
