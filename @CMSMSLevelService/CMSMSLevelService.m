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
            %   msms_cfg (CMSMSPepDeconvConfig)
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

            mapModification = readModifyInfo(cfg.mod_file_path);
            fixedModNameMass = obj.getModMassName(cfg.fixed_mod, mapModification);
            variableModNameMass = obj.getModMassName(cfg.variable_mod, mapModification);

            cMgfDatasetIO = CMgfDatasetIO(cfg.spec_dir_path);
            pepProtService = CPepProtService(cfg.fasta_file_path, cfg.regular_express, cfg.filtered_res_file_path);
            cMsFileMapper = CMsFileMapper(cfg.spec_dir_path); %#ok<NASGU>

            fin = fopen(cfg.pep_spec_file_path,'r');
            if 0>=fin
                error(['Can not open file: ',cfg.pep_spec_file_path]);
            end
            file_total_length = dir(cfg.pep_spec_file_path).bytes;
            print_progress = CPrintProgress(file_total_length);
            fprintf('Quantifying at PSM level...')
            warning_message = [];

            CPathResolver.ensureDir(cfg.output_dir_path);

            % TODO: Set the output path for MSMS level results
            each_PSM_results_path = CPathResolver.resolveFilePath(cfg.output_dir_path, 'report_msms.txt', '');
            msms_result = CMS2Result();
            may_fp_report_path = CPathResolver.resolveFilePath(cfg.output_dir_path, 'report_spectra_may_FP.txt', '');
            fo_may_FP = fopen(may_fp_report_path,'w');
            if fo_may_FP <= 0
                error(['Cannot open the the report file ', may_fp_report_path]);
            end
            strLine = fgetl(fin);
            str = regexp(strLine,'\t','split');
            pepSeq = str{1};
            if_wrote_peptide = false;

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
                    if_wrote_peptide = false;
                    pepSeq = str{1};
                else
                    % a spectrum for an old peptide
                    dataset_name = str{1};
                    spec_name = str{2};
                    [isProtN,isProtC] = pepProtService.getWhetherProtNC(pepSeq);
                    eachSpecPipeline = CMS2SpectrumPipeline(pepSeq,isProtN,isProtC, ...
                        cMgfDatasetIO,dataset_name,spec_name,fixedModNameMass, ...
                        variableModNameMass,cms2_cfg);

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

                        [bSuccess,cstrIMP,abundance,ionTypePosCharge,ionIntens,frageff,warning_msg,is_X_not_full_column_rank] = ...
                            eachSpecPipeline.processSpectrumWithContext(peptideCtx, spectrumCtx);
                    catch ME
                        bSuccess = false;
                        cstrIMP = [];
                        abundance = [];
                        ionTypePosCharge = [];
                        ionIntens = [];
                        frageff = [];
                        is_X_not_full_column_rank = false;
                        warning_msg = ['[CMS2] ', ME.identifier, ': ', ME.message, '\n'];
                    end

                    warning_message = [warning_message, warning_msg]; %#ok<AGROW>
                    if bSuccess
                        if is_X_not_full_column_rank
                            fprintf(fo_may_FP,'%s\t%s\n',dataset_name,spec_name);
                        end

                        if obj.m_is_record_fragment_information
                            obj.appendFragmentInformation(ionTypePosCharge, ionIntens, frageff);
                        end

                        if ~if_wrote_peptide
                            msms_result.addPeptide(pepSeq);
                            if_wrote_peptide = true;
                        end
                        msms_result.addSpectrum(dataset_name,spec_name);
                        imp_idx_nonzero = find(abundance~=0);
                        for idx = 1:length(imp_idx_nonzero)
                            msms_result.addPeptidoform(cstrIMP{imp_idx_nonzero(idx)}, ...
                                abundance(imp_idx_nonzero(idx)));
                        end
                    end
                end
            end

            CMS2ResultIO.write(msms_result, each_PSM_results_path);
            fclose(fo_may_FP);
            fclose(fin);
            print_progress.last_update();
            fprintf('done.\n');
            if warning_message
                fprintf(warning_message);
            end
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

        function [modNameMass]=getModMassName(~,modificationTypes,mapModification)
            modNameMass=[];
            if isempty(modificationTypes)
                return
            end
            S_modificationTypes=regexp(modificationTypes,';','split');
            modNameMass=cell(length(S_modificationTypes),3);
            for i=1:length(S_modificationTypes)
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
            modNameMass(cellfun(@isempty,modNameMass(:,1)),:)=[];
        end
    end
end
