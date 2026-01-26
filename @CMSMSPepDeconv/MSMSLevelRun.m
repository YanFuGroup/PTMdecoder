function obj = MSMSLevelRun(obj, is_record_fragment_information)
% MSMSLevelRun - Quantify the MSMS level results
% Input:
%   obj - The CMSMSPepDeconv object
%   is_record_fragment_information - Whether to record the fragment information
% Output:
%   obj - The CMSMSPepDeconv object

% Load protein sequences from fasta file
if isempty(obj.CPepProtService)
    obj.CPepProtService = CPepProtService(obj.m_fastaFile, obj.m_regular_express, obj.m_filtered_res_file_path);
end

% Indexing the mgf
if isempty(obj.m_cMgfDatasetIO)
    obj.m_cMgfDatasetIO = CMgfDatasetIO;
    obj.m_cMgfDatasetIO.Init(obj.m_specPath);
    obj.m_cMgfDatasetIO.SetMap();
    obj.m_cMgfDatasetIO.SetFidmap();
    need_release_mgf_index = true;
else
    need_release_mgf_index = false;
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

mkdir(obj.m_outputDir);

each_PSM_results_path = fullfile(obj.m_outputDir,'report_msms.txt');
fout = fopen(each_PSM_results_path,'w');
if fout <= 0
    error(['Cannot open the the report file ', each_PSM_results_path]);
end
fo_may_FP = fopen(fullfile(obj.m_outputDir,'report_spectra_may_FP.txt'),'w');
if fo_may_FP <= 0
    error(['Cannot open the the report file ', fullfile(obj.m_outputDir,'report_spectra_may_FP.txt')]);
end
strLine = fgetl(fin);
str = regexp(strLine,'\t','split');
pepSeq = str{1}; % record the current peptide sequence
if_wrote_peptide = false;
if_the_first = true;
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
        [isProtN,isProtC] = obj.CPepProtService.getWhetherProtNC(pepSeq);
        eachSpecIMSLQ = CEachSpectrumLocQuant(pepSeq,isProtN,isProtC, ...
            obj.m_cMgfDatasetIO,str{1},str{2},obj.m_fixedModNameMass, ...
            obj.m_variableModNameMass,obj.m_model,obj.m_method,obj.m_lambda, ...
            obj.m_ms1_tolerance,obj.m_ms2_tolerance,obj.m_alpha,obj.m_resFilterThres,...
            obj.m_ionTypes,obj.m_enzyme);
        [bSuccess,cstrPepIso,abundance,ionTypePosCharge,ionIntens,frageff, ...
            warning_msg,is_X_full_column_rank] = eachSpecIMSLQ.runEach();
        warning_message = [warning_message, warning_msg]; %#ok<AGROW> 
        if bSuccess
            if is_X_full_column_rank
                fprintf(fo_may_FP,'%s\t%s\n',str{1},str{2});
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
                if ~if_the_first
                    fprintf(fout,'\n\n');
                end
                fprintf(fout,'P\t%s\n',pepSeq);
                if_wrote_peptide = true;
                if_the_first = false;   
            end
            fprintf(fout,'S\t%s\t%s\n',str{1},str{2});
            idxNonZero = find(abundance~=0);
            for idx = 1:length(idxNonZero)
                fprintf(fout,'%s\t%f\n',cstrPepIso{idxNonZero(idx)},...
                    abundance(idxNonZero(idx)));
            end
        end
    end
end
fclose(fout);
fclose(fo_may_FP);
fclose(fin);
if need_release_mgf_index
    obj.m_cMgfDatasetIO.CloseAllFile();
end
print_progress.last_update();
fprintf('done.\n');
if warning_message
    fprintf(warning_message);
end
end