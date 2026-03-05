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
            %       - msms_res_path (1 x 1 char/string)
            %           input MSMS result path (default: cfg.msms_res_path or report_msms.txt)
            %       - peptide_quant_res_path (1 x 1 char/string)
            %           initial peptide quant result path (default: report_peptide_all.txt)
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

            CPathResolver.ensureDir(cfg.output_dir_path);

            msms_res_path = CPathResolver.resolveFilePath(cfg.output_dir_path, 'report_msms.txt', ...
                CStructOptionUtils.get(align_options, 'msms_res_path', cfg.msms_res_path));
            msms_result = CMS2ResultIO.read(msms_res_path);

            quant_output_path = CPathResolver.resolveFilePath(cfg.output_dir_path, 'report_peptide_all_primary.txt', ...
                CStructOptionUtils.get(align_options, 'peptide_quant_res_path', ''));

            rawIdentManagers = cell(1, length(msms_result.Peptides));
            for i_pep = 1:length(msms_result.Peptides)
                rawIdentManagers{i_pep} = obj.buildRawIdentManagerFromSpectrumList(msms_result.Peptides(i_pep).spectrum_list);
            end

            proc_cfg = CIMPProcessingExecutorConfig(struct( ...
                'ms12DatasetIO', obj.m_cMs12DatasetIO, ...
                'ms1_tolerance', cfg.ms1_tolerance, ...
                'minMSMSnum', cfg.min_MSMS_num, ...
                'alpha', cfg.alpha, ...
                'resFilterThres', 0));  % cfg.result_filter_threshold
            proc_executor = CIMPProcessingExecutor(proc_cfg);

            quant_report = CIMPQuantReport();
            print_progress = CPrintProgress(length(msms_result.Peptides));
            fprintf('Quantifying at peptide level before alignment...')
            for i_pep = 1:length(msms_result.Peptides)
                print_progress = print_progress.update_show(i_pep);
                peptide_sequence = msms_result.Peptides(i_pep).peptide_sequence;
                cell_prot_name_pos = obj.m_pepProtService.get_protein_name_pos(peptide_sequence);
                rawIdentManager = rawIdentManagers{i_pep};
                block = proc_executor.quantifyPeptideBlock(cell_prot_name_pos, rawIdentManager);
                quant_report = quant_report.append_block(block);
            end
            print_progress.last_update();
            fprintf('done.\n');

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
                cfg.filtered_res_file_path, rawIdentManagers, base_pep_rtrange_map);

            align_executor.writeAlignmentReport(align_report, cfg.alignment_report_path);

            output_path = cfg.requant_output_path;
            report = CIMPQuantReport();
            print_progress = CPrintProgress(length(msms_result.Peptides));

            fprintf('Re-quantifying at peptide level (aligned)...')
            for i_pep = 1:length(msms_result.Peptides)
                print_progress = print_progress.update_show(i_pep);
                peptide_sequence = msms_result.Peptides(i_pep).peptide_sequence;
                cell_prot_name_pos = obj.m_pepProtService.get_protein_name_pos(peptide_sequence);
                rawIdentManager = rawIdentManagers{i_pep};
                block = proc_executor.requantifyPeptideBlock(cell_prot_name_pos, rawIdentManager, pep_rtrange_map);
                report = report.append_block(block);
            end
            print_progress.last_update();
            fprintf('done.\n');

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

    end
end
