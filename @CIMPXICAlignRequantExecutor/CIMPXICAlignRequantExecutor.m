classdef CIMPXICAlignRequantExecutor < handle
    % Execution layer for quant -> XIC align -> requant workflow

    properties (Access=private)
        m_ms12DatasetIO
        m_ms1_tolerance
        m_minMSMSnum
        m_resFilterThres
        m_aligner
        m_align_strategy
        m_align_options

        m_pair_models
    end

    methods
        function obj = CIMPXICAlignRequantExecutor(ms12DatasetIO, ms1_tolerance, minMSMSnum, ...
            resFilterThres, aligner, align_strategy, align_options)
            % Construct a quant-align-requant executor.
            % Input:
            %   ms12DatasetIO (CMS12DatasetIO)
            %       MS1/MS2 dataset IO
            %   ms1_tolerance (struct)
            %       MS1 tolerance settings
            %   minMSMSnum (1 x 1 double/int)
            %       Minimum MSMS count for XIC preprocessing
            %   resFilterThres (1 x 1 double)
            %       Ratio filter threshold
            %   aligner (CXICAligner)
            %       Aligner instance
            %   align_strategy (CRunAlignStrategy)
            %       Alignment strategy
            %   align_options (struct, optional)
            %       Alignment options (paths and parameters):
            %         - min_psm (1 x 1 double)
            %             Minimum PSM count per peptide to form anchors. Default: 3.
            %         - num_bins (1 x 1 double)
            %             Number of bins for local alignment offsets. Default: 5.
            %         - min_per_bin (1 x 1 double)
            %             Minimum anchors per bin for local offset. Default: 5.
            %         - outlier_k (1 x 1 double)
            %             Outlier removal k (MAD/STD) for transform fitting. Default: 3.
            %         - outlier_method (char/string)
            %             Outlier method ('mad' or 'std'). Default: 'mad'.
            %         - rt_sigma (1 x 1 double)
            %             RT Gaussian sigma (minutes) for peak selection. Default: 0.5.
            %         - max_rt_residual (1 x 1 double)
            %             Max allowed RT residual (minutes) for peak pairing. Default: model.rt_residual_threshold.
            %         - dead_time_min (1 x 1 double)
            %             Min allowed RT start (minutes) for peaks. Default: [] (disabled).
            % Output:
            %   obj (CIMPXICAlignRequantExecutor)
            %       Executor instance
            if nargin == 1 && isa(ms12DatasetIO, 'CIMPXICAlignRequantExecutorConfig')
                cfg = ms12DatasetIO;
                obj.m_ms12DatasetIO = cfg.ms12DatasetIO;
                obj.m_ms1_tolerance = cfg.ms1_tolerance;
                obj.m_minMSMSnum = cfg.minMSMSnum;
                obj.m_resFilterThres = cfg.resFilterThres;
                obj.m_aligner = cfg.aligner;
                obj.m_align_strategy = cfg.align_strategy;
                obj.m_align_options = cfg.align_options;
                return;
            end

            obj.m_ms12DatasetIO = ms12DatasetIO;
            obj.m_ms1_tolerance = ms1_tolerance;
            obj.m_minMSMSnum = minMSMSnum;
            obj.m_resFilterThres = resFilterThres;
            obj.m_aligner = aligner;
            obj.m_align_strategy = align_strategy;
            if nargin < 8
                align_options = struct();
            end
            obj.m_align_options = align_options;
        end
    end

    methods (Access=public)
        [pep_rtrange_map, report] = buildAlignedRtRangeMap(obj, fdr_filtered_result_path, rawIdentManagers, base_pep_rtrange_map, base_groups_by_peptide)

        writeAlignmentReport(obj, report, output_path)
    end

    methods (Access=private)
        raw_imps_map = buildRawImpsMap(obj, rawIdentManager)

        raw_imps_map = buildRawImpsMapFromBaseGroups(obj, base_groups)

        raw_names = getRawNamesFromRawIdentManagers(obj, rawIdentManagers)
    end
end
