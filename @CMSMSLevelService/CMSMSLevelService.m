classdef CMSMSLevelService < handle
    % Stage service for MSMS-level quantification.
    
    properties (Access = private)
        m_msms_cfg
        m_is_record_fragment_information
        m_matFragInfo
        m_matFragEff
        m_matFragIntens
    end
    
    methods
        function obj = CMSMSLevelService(msms_cfg)
            % Input:
            %   msms_cfg (struct)
            %       config for MSMS processor
            if nargin < 1 || isempty(msms_cfg)
                error('CMSMSLevelService:MissingConfig', ...
                    'msms_cfg must be provided.');
            end
            obj.m_msms_cfg = msms_cfg;
            obj.m_is_record_fragment_information = false;
            obj.m_matFragInfo = [];
            obj.m_matFragEff = [];
            obj.m_matFragIntens = [];
        end
        
        function setRecordFragmentInformation(obj, enabled)
            % Enable/disable fragment-information recording.
            % Input:
            %   enabled (1 x 1 logical)
            %       true to enable recording; false to disable
            if nargin < 2
                enabled = true;
            end
            obj.m_is_record_fragment_information = logical(enabled);
        end
        
        function frag = getFragmentInformation(obj)
            % Get recorded fragment-information buffers.
            % Output:
            %   frag (struct)
            %       - matFragInfo
            %       - matFragEff
            %       - matFragIntens
            frag = struct( ...
                'matFragInfo', obj.m_matFragInfo, ...
                'matFragEff', obj.m_matFragEff, ...
                'matFragIntens', obj.m_matFragIntens);
        end
        
        function run(obj)
            % Run MSMS-level stage.
            cfg = obj.m_msms_cfg;
            
            [fixedModNameMass, variableModNameMass] = CModificationRegistry.fromConfig(cfg);
            
            cMgfDatasetIO = CMgfDatasetIO(cfg.spec_dir_path);
            pepProtService = CPepProtService(cfg.fasta_file_path, cfg.regular_express, cfg.filtered_res_file_path);
            cMsFileMapper = CMsFileMapper(cfg.spec_dir_path); %#ok<NASGU>
            
            fin = fopen(cfg.pep_spec_file_path,'r');
            if 0>=fin
                error(['Can not open file: ',cfg.pep_spec_file_path]);
            end
            file_total_length = dir(cfg.pep_spec_file_path).bytes;
            print_progress = CPrintProgress(file_total_length, 'MSMS level processing');
            CLogger.info('Quantifying at PSM level...');
            
            CPathResolver.ensureDir(cfg.output_dir_path);
            
            % TODO: Set the output path for MSMS level results
            each_PSM_results_path = CPathResolver.resolveFilePath(cfg.output_dir_path, 'report_msms.txt', '');
            may_fp_report_path = CPathResolver.resolveFilePath(cfg.output_dir_path, 'report_spectra_may_FP.txt', '');
            fo_may_FP = fopen(may_fp_report_path,'w');
            if fo_may_FP <= 0
                error(['Cannot open the the report file ', may_fp_report_path]);
            end
            strLine = fgetl(fin);
            str = regexp(strLine,'\t','split');
            pepSeq = str{1};
            stage1_records = struct( ...
                'dataset_name', {}, ...
                'spec_name', {}, ...
                'pepSeq', {}, ...
                'bSuccess', {}, ...
                'isShortcut', {}, ...
                'is_X_not_full_column_rank', {}, ...
                'abundance', {}, ...
                'cstrIMP', {}, ...
                'reported_imp_write_indices', {}, ...
                'solver_diag', {}, ...
                'noise_model_fit_inputs', {}, ...
                'stability_cache', {});
            stage1_total_count = 0;
            stage1_success_count = 0;
            stage1_shortcut_count = 0;
            
            % Build and reuse static CMS2 pipeline config for all spectra in this run
            cms2_cfg = CMS2SpectrumPipelineConfig(struct( ...
                'model', cfg.model, ...
                'method', cfg.method, ...
                'lambda', cfg.lambda, ...
                'ms1_tolerance', cfg.ms1_tolerance, ...
                'ms2_tolerance', cfg.ms2_tolerance, ...
                'alpha', cfg.alpha, ...
                'resFilterThres', cfg.result_filter_threshold, ...
                'ionTypes', cfg.ion_types, ...
                'enzyme', struct('name', cfg.enzyme_name, 'limits', cfg.enzyme_limits), ...
                'case_penalty_intens', cfg.case_penalty_intens, ...
                'grid_penalty_intens', cfg.grid_penalty_intens, ...
                'case_OLS_intens_weight', cfg.case_OLS_intens_weight));
            eachSpecPipeline = CMS2SpectrumPipeline(fixedModNameMass, variableModNameMass, cms2_cfg);
            
            while ~feof(fin)
                strLine = fgetl(fin);
                now_bytes = ftell(fin);
                print_progress = print_progress.update_show(now_bytes);
                if isempty(strtrim(strLine))
                    continue;
                end
                str = regexp(strLine,'\t','split');
                if length(str)==1 || isempty(str{2})
                    % a new peptide
                    pepSeq = str{1};
                else
                    % a spectrum for an old peptide
                    dataset_name = str{1};
                    spec_name = str{2};
                    [isProtN,isProtC] = pepProtService.getWhetherProtNC(pepSeq);
                    
                    try
                        [expPeaks,iCharge,precursorMZ] = cMgfDatasetIO.read_oneSpec(dataset_name,spec_name);
                        
                        peptideCtx = struct( ...
                            'pepSeq', pepSeq, ...
                            'isProtN', isProtN, ...
                            'isProtC', isProtC, ...
                            'fixedModNameMass', {fixedModNameMass}, ...
                            'variableModNameMass', {variableModNameMass});
                        
                        spectrumCtx = struct( ...
                            'datasetName', dataset_name, ...
                            'specName', spec_name, ...
                            'expPeaks', expPeaks, ...
                            'iCharge', iCharge, ...
                            'precursorMZ', precursorMZ);
                        
                        [bSuccess,cstrIMP,abundance,ionTypePosCharge,ionIntens,frageff, ...
                            is_X_not_full_column_rank,solver_diag,noise_model_fit_inputs,stability_cache] = ...
                            eachSpecPipeline.runBaselineSpectrumStage(peptideCtx, spectrumCtx);
                    catch ME
                        bSuccess = false;
                        cstrIMP = cell(0,1);
                        abundance = [];
                        ionTypePosCharge = [];
                        ionIntens = [];
                        frageff = [];
                        is_X_not_full_column_rank = false;
                        solver_diag = struct( ...
                            'is_X_not_full_column_rank', false, ...
                            'jaccard_stability', NaN, ...
                            'support_frequency', [], ...
                            'abundance_mad', [], ...
                            'reported_imp_indices', [], ...
                            'num_successful_resamples', 0);
                        noise_model_fit_inputs = struct( ...
                            'filteredOutExpPeakCount', 0, ...
                            'filteredOutExpPeakSqSum', 0, ...
                            'matchedExpPeaks', zeros(0, 3), ...
                            'fittedMatchedPeakIntensities', zeros(0, 1));
                        stability_cache = struct( ...
                            'vNonRedunTheoryIonMz', [], ...
                            'matchedExpPeaks', zeros(0, 3), ...
                            'massArrangement', [], ...
                            'abundance', [], ...
                            'fittedMatchedPeakIntensities', zeros(0, 1), ...
                            'cstrIMP', {{}}, ...
                            'solver_diag', solver_diag);
                        CLogger.debug('[CMSMSLevelService:run] Failed spectrum (%s/%s) - %s: %s', ...
                            dataset_name, spec_name, ME.identifier, ME.message);
                    end
                    
                    stage1_total_count = stage1_total_count + 1;
                    stage1_is_shortcut = bSuccess && (isempty(stability_cache.matchedExpPeaks) || size(stability_cache.massArrangement,1) <= 1);
                    stage1_record = struct( ...
                        'dataset_name', dataset_name, ...
                        'spec_name', spec_name, ...
                        'pepSeq', pepSeq, ...
                        'bSuccess', bSuccess, ...
                        'isShortcut', stage1_is_shortcut, ...
                        'is_X_not_full_column_rank', is_X_not_full_column_rank, ...
                        'abundance', abundance, ...
                        'cstrIMP', {cstrIMP}, ...
                        'reported_imp_write_indices', [], ...
                        'solver_diag', solver_diag, ...
                        'noise_model_fit_inputs', noise_model_fit_inputs, ...
                        'stability_cache', stability_cache);
                    if bSuccess
                        stage1_success_count = stage1_success_count + 1;
                        if stage1_is_shortcut
                            stage1_shortcut_count = stage1_shortcut_count + 1;
                        end
                        
                        if is_X_not_full_column_rank
                            fprintf(fo_may_FP,'%s\t%s\n',dataset_name,spec_name);
                        end
                        
                        if obj.m_is_record_fragment_information
                            obj.appendFragmentInformation(ionTypePosCharge, ionIntens, frageff);
                        end
                        
                        imp_idx_nonzero = find(abundance~=0);
                        stage1_record.reported_imp_write_indices = imp_idx_nonzero;
                    end
                    stage1_records(end + 1) = stage1_record;
                end
            end
            
            CLogger.info('[CMSMSLevelService:run] Stage-1 baseline collection done. total=%d, success=%d, shortcut=%d.', ...
                stage1_total_count, stage1_success_count, stage1_shortcut_count);
            
            stability_options = struct();
            has_stability_options = false;
            if isfield(cfg, 'stability_options') && isstruct(cfg.stability_options)
                stability_options = cfg.stability_options;
                has_stability_options = true;
            end
            
            noise_model = [];
            stage2_total_candidates = 0;
            stage2_success_count = 0;
            stage2_failed_count = 0;
            
            if has_stability_options
                if isempty(stage1_records)
                    eligible_noise_fit_records = struct([]);
                else
                    eligible_mask = arrayfun(@(rec) rec.bSuccess && ~rec.isShortcut && ...
                        isfield(rec.noise_model_fit_inputs, 'matchedExpPeaks') && ...
                        ~isempty(rec.noise_model_fit_inputs.matchedExpPeaks), stage1_records);
                    eligible_noise_fit_records = [stage1_records(eligible_mask).noise_model_fit_inputs];
                end
                
                if isempty(eligible_noise_fit_records)
                    CLogger.warn('[CMSMSLevelService:run] No eligible spectra for dataset-level noise fitting. Stage-2 stability is skipped.');
                else
                    try
                        noise_model = CMS2QuantSolver.estimateDatasetNoiseModel(eligible_noise_fit_records);
                    catch ME
                        noise_model = [];
                        CLogger.warn('[CMSMSLevelService:run] estimateDatasetNoiseModel failed: [%s] %s. Stage-2 stability is skipped.', ...
                            ME.identifier, ME.message);
                    end
                end
            else
                CLogger.info(['[CMSMSLevelService:run] Stage-2 stability is disabled because cfg.stability_options is not configured. ', ...
                    'Set n_resamples in the parameter file to enable Stage-2 stability.']);
            end
            
            if ~isempty(noise_model) && has_stability_options
                solver_cfg = struct( ...
                    'model', cfg.model, ...
                    'method', cfg.method, ...
                    'lambda', cfg.lambda, ...
                    'case_penalty_intens', cfg.case_penalty_intens, ...
                    'grid_penalty_intens', cfg.grid_penalty_intens, ...
                    'case_OLS_intens_weight', cfg.case_OLS_intens_weight);
                
                for idxRecord = 1:numel(stage1_records)
                    rec = stage1_records(idxRecord);
                    if ~(rec.bSuccess && ~rec.isShortcut)
                        continue;
                    end
                    
                    stage2_total_candidates = stage2_total_candidates + 1;
                    
                    cache = rec.stability_cache;
                    if isempty(cache.vNonRedunTheoryIonMz) || isempty(cache.matchedExpPeaks) || isempty(cache.massArrangement) || isempty(cache.abundance)
                        stage2_failed_count = stage2_failed_count + 1;
                        CLogger.warn('[CMSMSLevelService:run] Stage-2 skipped for %s/%s due to incomplete stability cache.', ...
                            rec.dataset_name, rec.spec_name);
                        continue;
                    end
                    
                    try
                        stability_diag = CMS2QuantSolver.estimateStability( ...
                            cache.vNonRedunTheoryIonMz, ...
                            cache.matchedExpPeaks, ...
                            cache.massArrangement, ...
                            solver_cfg, ...
                            cache.abundance, ...
                            cache.fittedMatchedPeakIntensities, ...
                            noise_model, ...
                            stability_options);
                        
                        rec.solver_diag.jaccard_stability = stability_diag.jaccard_stability;
                        rec.solver_diag.support_frequency = stability_diag.support_frequency;
                        rec.solver_diag.abundance_mad = stability_diag.abundance_mad;
                        rec.solver_diag.reported_imp_indices = stability_diag.reported_imp_indices;
                        rec.solver_diag.num_successful_resamples = stability_diag.num_successful_resamples;
                        rec.stability_cache.solver_diag = rec.solver_diag;
                        
                        stage1_records(idxRecord) = rec;
                        stage2_success_count = stage2_success_count + 1;
                    catch ME
                        stage2_failed_count = stage2_failed_count + 1;
                        CLogger.warn('[CMSMSLevelService:run] estimateStability failed for %s/%s: [%s] %s', ...
                            rec.dataset_name, rec.spec_name, ME.identifier, ME.message);
                    end
                end

                CLogger.info('[CMSMSLevelService:run] Stage-2 stability done. candidates=%d, success=%d, failed=%d.', ...
                stage2_total_candidates, stage2_success_count, stage2_failed_count);
            end
            
            msms_result = CMS2Result();
            for idxRecord = 1:numel(stage1_records)
                rec = stage1_records(idxRecord);
                if ~rec.bSuccess
                    continue;
                end
                
                msms_result.addOrSelectPeptide(rec.pepSeq);
                
                spectrumMeta = struct('jaccard_stability', rec.solver_diag.jaccard_stability);
                msms_result.addSpectrum(rec.dataset_name, rec.spec_name, spectrumMeta);
                
                support_full = NaN(size(rec.abundance));
                mad_full = NaN(size(rec.abundance));
                reported_imp_indices = rec.solver_diag.reported_imp_indices;
                support_frequency = rec.solver_diag.support_frequency;
                abundance_mad = rec.solver_diag.abundance_mad;
                if ~isempty(reported_imp_indices) && ~isempty(support_frequency)
                    if numel(reported_imp_indices) ~= numel(support_frequency)
                        CLogger.error(['[CMSMSLevelService:InconsistentSupportFrequencyLength] ', ...
                            'reported_imp_indices length (%d) must equal support_frequency length (%d) for %s/%s.'], ...
                            numel(reported_imp_indices), numel(support_frequency), rec.dataset_name, rec.spec_name);
                    end
                    out_of_range_mask = reported_imp_indices < 1 | reported_imp_indices > numel(support_full);
                    if any(out_of_range_mask)
                        CLogger.error(['[CMSMSLevelService:InvalidReportedImpIndex] ', ...
                            'reported_imp_indices contains out-of-range entries for %s/%s.'], ...
                            rec.dataset_name, rec.spec_name);
                    end
                    support_full(reported_imp_indices) = support_frequency;
                end
                if ~isempty(reported_imp_indices) && ~isempty(abundance_mad)
                    if numel(reported_imp_indices) ~= numel(abundance_mad)
                        CLogger.error(['[CMSMSLevelService:InconsistentAbundanceMadLength] ', ...
                            'reported_imp_indices length (%d) must equal abundance_mad length (%d) for %s/%s.'], ...
                            numel(reported_imp_indices), numel(abundance_mad), rec.dataset_name, rec.spec_name);
                    end
                    out_of_range_mask = reported_imp_indices < 1 | reported_imp_indices > numel(mad_full);
                    if any(out_of_range_mask)
                        CLogger.error(['[CMSMSLevelService:InvalidReportedImpIndex] ', ...
                            'reported_imp_indices contains out-of-range entries for %s/%s.'], ...
                            rec.dataset_name, rec.spec_name);
                    end
                    mad_full(reported_imp_indices) = abundance_mad;
                end
                
                imp_idx_nonzero = rec.reported_imp_write_indices;
                for idxImp = 1:length(imp_idx_nonzero)
                    imp_idx = imp_idx_nonzero(idxImp);
                    peptidoformMeta = struct( ...
                        'support_frequency', support_full(imp_idx), ...
                        'abundance_mad', mad_full(imp_idx));
                    msms_result.addPeptidoform(rec.cstrIMP{imp_idx}, rec.abundance(imp_idx), peptidoformMeta);
                end
            end
            
            CMS2ResultIO.write(msms_result, each_PSM_results_path);
            fclose(fo_may_FP);
            fclose(fin);
            print_progress.last_update();
            CLogger.info('MSMS-level quantification done.');
        end
    end
    
    methods (Access = private)
        function appendFragmentInformation(obj, ionTypePosCharge, ionIntens, frageff)
            if isempty(frageff)
                % Skip the fragment information of this spectrum if it is empty
                %   (only one possible peptidoform and not have been deconvoluted)
                return;
            end
            
            if isempty(obj.m_matFragInfo)
                % Allocate space and add a column of data
                obj.m_matFragInfo = ionTypePosCharge;
                obj.m_matFragIntens = ionIntens;
                obj.m_matFragEff = frageff;
            else
                obj.m_matFragEff = [obj.m_matFragEff, zeros(size(obj.m_matFragEff,1),1)];
                obj.m_matFragIntens = [obj.m_matFragIntens, zeros(size(obj.m_matFragIntens,1),1)];
                for idxFrag = 1:size(ionTypePosCharge,1)
                    [bIsExist,idxMatFE] = ismember(ionTypePosCharge(idxFrag,:), obj.m_matFragInfo, 'rows');
                    % Check whether this ion type exists, if it exists, record it directly, if not, allocate space and then record it
                    if bIsExist
                        obj.m_matFragEff(idxMatFE,end) = frageff(idxFrag);
                        obj.m_matFragIntens(idxMatFE,end) = ionIntens(idxFrag);
                    else
                        obj.m_matFragInfo = [obj.m_matFragInfo; ionTypePosCharge(idxFrag,:)];
                        obj.m_matFragIntens = [obj.m_matFragIntens; zeros(1,size(obj.m_matFragIntens,2))];
                        obj.m_matFragIntens(end,end) = ionIntens(idxFrag);
                        obj.m_matFragEff = [obj.m_matFragEff; zeros(1,size(obj.m_matFragEff,2))];
                        obj.m_matFragEff(end,end) = frageff(idxFrag);
                    end
                end
            end
        end
        
    end
end
