function obj = MSMSLevelRun(obj, is_record_fragment_information)
% MSMSLevelRun - Quantify the MSMS level results
% Input:
%   obj (CMSMSPepDeconv)
%       Processor instance
%   is_record_fragment_information (1 x 1 logical)
%       Whether to record fragment information
% Output:
%   obj (CMSMSPepDeconv)
%       Updated instance

% Load shared services lazily
[obj, ~] = obj.ensurePepProtService();

[obj, ~] = obj.ensureMsFileMapper();

[obj, need_release_mgf_index] = obj.ensureMgfDatasetIO();
if need_release_mgf_index
    cleanup_mgf = onCleanup(@() obj.m_cMgfDatasetIO.CloseAllFile()); %#ok<NASGU>
end


% Read the "peptide-spectrum" list and start analysis centered on peptides
fin = fopen(obj.m_pepSpecFile,'r');
if 0>=fin
    error(['Can not open file: ',obj.m_pepSpecFile]);
end
file_total_length = dir(obj.m_pepSpecFile).bytes;
print_progress = CPrintProgress(file_total_length);
fprintf('Quantifying at PSM level...')
warning_message = [];

if ~isfolder(obj.m_outputDir)
    mkdir(obj.m_outputDir);
end

each_PSM_results_path = fullfile(obj.m_outputDir,'report_msms.txt');
msms_result = CMS2Result();
fo_may_FP = fopen(fullfile(obj.m_outputDir,'report_spectra_may_FP.txt'),'w');
if fo_may_FP <= 0
    error(['Cannot open the the report file ', fullfile(obj.m_outputDir,'report_spectra_may_FP.txt')]);
end
strLine = fgetl(fin);
str = regexp(strLine,'\t','split');
pepSeq = str{1}; % record the current peptide sequence
if_wrote_peptide = false;

% Build and reuse static CMS2 pipeline config for all spectra in this run
cms2_cfg = CMS2SpectrumPipelineConfig(struct( ...
    'model', obj.m_model, ...
    'method', obj.m_method, ...
    'lambda', obj.m_lambda, ...
    'ms1_tolerance', obj.m_ms1_tolerance, ...
    'ms2_tolerance', obj.m_ms2_tolerance, ...
    'alpha', obj.m_alpha, ...
    'resFilterThres', obj.m_resFilterThres, ...
    'ionTypes', obj.m_ionTypes, ...
    'enzyme', obj.m_enzyme, ...
    'case_penalty_intens', obj.m_case_penalty_intens, ...
    'grid_penalty_intens', obj.m_grid_penalty_intens, ...
    'case_OLS_intens_weight', obj.m_case_OLS_intens_weight));

while ~feof(fin)
    strLine = fgetl(fin);
    now_bytes = ftell(fin);
    print_progress = print_progress.update_show(now_bytes);
    if isempty(strtrim(strLine))
        continue;
    end
    str = regexp(strLine,'\t','split');
    if length(str)==1 || isempty(str{2})
        % meet a new peptide
        if_wrote_peptide = false;
        pepSeq = str{1};
    else
        % meet a spectrum for an old peptide
        str = regexp(strLine,'\t','split');
        dataset_name = str{1};
        spec_name = str{2};
        [isProtN,isProtC] = obj.CPepProtService.getWhetherProtNC(pepSeq);
        eachSpecPipeline = CMS2SpectrumPipeline(pepSeq,isProtN,isProtC, ...
            obj.m_cMgfDatasetIO,dataset_name,spec_name,obj.m_fixedModNameMass, ...
            obj.m_variableModNameMass,cms2_cfg);

        try
            [expPeaks,iCharge,precursorMZ] = obj.m_cMgfDatasetIO.read_oneSpec( ...
                dataset_name,spec_name);

            peptideCtx = struct( ...
                'pepSeq', pepSeq, ...
                'isProtN', isProtN, ...
                'isProtC', isProtC, ...
                'fixedModNameMass', {obj.m_fixedModNameMass}, ...
                'variableModNameMass', {obj.m_variableModNameMass});

            spectrumCtx = struct( ...
                'datasetName', dataset_name, ...
                'specName', spec_name, ...
                'expPeaks', expPeaks, ...
                'iCharge', iCharge, ...
                'precursorMZ', precursorMZ);

            [bSuccess,cstrIMP,abundance,ionTypePosCharge,ionIntens,frageff, ...
                warning_msg,is_X_not_full_column_rank] = ...
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
            if is_record_fragment_information
                if isempty(frageff)
                    % Skip the fragment information of this spectrum if it is empty
                    %   (only one possible peptidoform and not have been deconvoluted)
                else
                    if isempty(obj.m_matFragInfo)% Allocate space and add a column of data
                        obj.m_matFragInfo=ionTypePosCharge;
                        obj.m_matFragIntens=ionIntens;
                        obj.m_matFragEff=frageff;
                    else
                        obj.m_matFragEff=[obj.m_matFragEff,zeros(size(obj.m_matFragEff,1),1)];
                        obj.m_matFragIntens=[obj.m_matFragIntens,zeros(size(obj.m_matFragIntens,1),1)];
                        for idxFrag=1:size(ionTypePosCharge,1)
                            [bIsExist,idxMatFE]=ismember(ionTypePosCharge(idxFrag,:), ...
                                obj.m_matFragInfo,'rows');
                            % Check whether this ion type exists, if it exists, record it directly, if not, allocate space and then record it
                            if bIsExist
                                obj.m_matFragEff(idxMatFE,end)=frageff(idxFrag);
                                obj.m_matFragIntens(idxMatFE,end)=ionIntens(idxFrag);
                            else
                                obj.m_matFragInfo=[obj.m_matFragInfo;ionTypePosCharge(idxFrag,:)];
                                obj.m_matFragIntens=[obj.m_matFragIntens;zeros(1,size(obj.m_matFragIntens,2))];
                                obj.m_matFragIntens(end,end) = ionIntens(idxFrag);
                                obj.m_matFragEff = [obj.m_matFragEff;zeros(1,size(obj.m_matFragEff,2))];
                                obj.m_matFragEff(end,end) = frageff(idxFrag);
                            end
                        end
                    end
                end
            end
            
            if ~if_wrote_peptide
                msms_result.addPeptide(pepSeq);
                if_wrote_peptide = true;
            end
            msms_result.addSpectrum(dataset_name,spec_name);
            imp_idx_nonzero = find(abundance~=0);
            for idx = 1:length(imp_idx_nonzero)
                msms_result.addPeptidoform(cstrIMP{imp_idx_nonzero(idx)},...
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