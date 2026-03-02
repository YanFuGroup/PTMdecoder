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
            print_progress = CPrintProgress(length(msms_result.Peptides));

            pipeline_cfg = CIMPProcessingExecutorConfig(struct( ...
                'ms12DatasetIO', cMs12DatasetIO, ...
                'ms1_tolerance', cfg.ms1_tolerance, ...
                'minMSMSnum', cfg.min_MSMS_num, ...
                'alpha', cfg.alpha, ...
                'resFilterThres', cfg.result_filter_threshold));
            executor = CIMPProcessingExecutor(pipeline_cfg);
            pipeline = CPeptideLevelPipeline(executor);
            stats_cleanup = onCleanup(@() CIMPQuantStats.rt_sorted_stats('flush', ...
                CPathResolver.resolveFilePath(cfg.output_dir_path, 'rt_sorted_stats.mat', '')));

            deps = struct( ...
                'getProfilesFunc', @(dataset_name, spectrum_name) obj.getProfiles( ...
                    cMgfDatasetIO, cMs12DatasetIO, cMsFileMapper, cfg.ms1_tolerance, cfg.spec_dir_path, dataset_name, spectrum_name), ...
                'fixedModNameMass', {fixedModNameMass}, ...
                'variableModNameMass', {variableModNameMass});

            fprintf('Quantifying at peptide level...')
            for idx_psf = 1:length(msms_result.Peptides)
                print_progress = print_progress.update_show(idx_psf);
                peptide_sequence = msms_result.Peptides(idx_psf).peptide_sequence;
                cell_prot_name_pos = pepProtService.get_protein_name_pos(peptide_sequence);
                rawIdentManager = CPeptideRawIdentAssembler.buildFromSpectrumList( ...
                    msms_result.Peptides(idx_psf).spectrum_list, deps);
                block = pipeline.quantifyPeptideBlock(cell_prot_name_pos, rawIdentManager);
                report = report.append_block(block);
            end
            print_progress.last_update();
            fprintf('done.\n');

            CIMPQuantResultIO.write(report, each_peptide_results_path);
        end
    end

    methods (Access = private)
        function [isorts,c_ref_isointens,cur_mz,cur_ch] = getProfiles(~, cMgfDatasetIO, cMs12DatasetIO, cMsFileMapper, ms1_tolerance, specPath, mgf_name, spectrum_name)
            % Read MS1 profile around one MS2 spectrum precursor.
            % Input:
            %   cMgfDatasetIO (CMgfDatasetIO)
            %       MGF dataset reader
            %   cMs12DatasetIO (CMS12DatasetIO)
            %       MS1/MS2 dataset reader
            %   cMsFileMapper (CMsFileMapper)
            %       MGF-to-MS1 mapping helper
            %   ms1_tolerance (struct)
            %       MS1 tolerance (fields: value, isppm)
            %   specPath (1 x 1 char/string)
            %       spectrum directory path
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
            [~, cur_ch, cur_mz] = cMgfDatasetIO.read_oneSpec(mgf_name,spectrum_name);

            mgf_stem = erase(mgf_name,'.mgf');
            if isempty(cMsFileMapper)
                mapper = CMsFileMapper(specPath);
            else
                mapper = cMsFileMapper;
            end
            ms12_stem = mapper.get_ms1_stem(mgf_stem);

            MS2_index = cMs12DatasetIO.m_mapNameMS2Index(ms12_stem);
            idx_cur_scan = MS2_index(:,2)==MS2ScanI;
            MS1Scan = MS2_index(idx_cur_scan,1);
            MS1_index = cMs12DatasetIO.m_mapNameMS1Index(ms12_stem);
            MS1_peaks = cMs12DatasetIO.m_mapNameMS1Peaks(ms12_stem);
            index_starts_MS1 = [1;MS1_index(1:size(MS1_index,1),3)];
            ino = find(MS1_index(:,1)==MS1Scan);
            isorts = MS1_index(ino,2);
            IX = index_starts_MS1(ino):index_starts_MS1(ino+1)-1;
            mz = MS1_peaks(IX,1);
            inten = MS1_peaks(IX,2);
            if ms1_tolerance.isppm
                ptol = ms1_tolerance.value*cur_mz*1e-6;
            else
                ptol = ms1_tolerance.value;
            end
            c_ptol = min([ptol,0.3]);
            c_ref_isointens = max(inten(abs(mz-cur_mz)<c_ptol));
            if isempty(c_ref_isointens)
                c_ref_isointens = 0;
            end
        end
    end
end
