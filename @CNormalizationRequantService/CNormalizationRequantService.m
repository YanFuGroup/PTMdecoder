classdef CNormalizationRequantService < handle
    % Stage service for normalization-peptide re-quantification.

    properties (Access = private)
        m_cfg
    end

    methods
        function obj = CNormalizationRequantService(cfg)
            % Input:
            %   cfg (struct)
            %       - msms_cfg (struct, minimal normalization config)
            %           required fields:
            %           spec_dir_path, ms1_tolerance, alpha,
            %           result_filter_threshold, output_dir_path,
            %           min_MSMS_num, checked_peptides_res_path
            %       - output_file_name (1 x 1 char/string)
            if nargin < 1 || isempty(cfg)
                error('CNormalizationRequantService:MissingConfig', ...
                    'cfg must be provided.');
            end
            obj.m_cfg = cfg;
        end

        function run(obj)
            % Run normalization peptide re-quantification stage.
            cfg = obj.m_cfg;
            msms_cfg = cfg.msms_cfg;

            checked_pep_path = CPathResolver.resolveFilePath(msms_cfg.output_dir_path, ...
                'peptide4normalization_checked.txt', msms_cfg.checked_peptides_res_path);

            report = CIMPQuantResultIO.read(checked_pep_path);
            pep_rtrange_map = report.build_pep_rtrange_map();
            pep_prot_map = report.build_pep_prot_map();
            if isempty(pep_rtrange_map)
                error(['The checked normalization peptide result file "', checked_pep_path, '" is empty!']);
            end

            CPathResolver.ensureDir(msms_cfg.output_dir_path);

            % TODO: Set the output path for requantified peptide level results
            output_path = CPathResolver.resolveFilePath(msms_cfg.output_dir_path, cfg.output_file_name, '');
            report_requant = CIMPQuantReport();

            cMs12DatasetIO = CMS12DatasetIO(msms_cfg.spec_dir_path, msms_cfg.ms1_tolerance);

            total_records = 0;
            for idx_block = 1:numel(report.blocks)
                total_records = total_records + numel(report.blocks(idx_block).records);
            end
            if total_records == 0
                CLogger.warn('The file "%s" is empty', checked_pep_path);
            end
            print_progress = CPrintProgress(max(total_records, 1), 'normalization_peptide_requantification');
            pipeline_cfg = CIMPProcessingExecutorConfig(struct( ...
                'ms12DatasetIO', cMs12DatasetIO, ...
                'ms1_tolerance', msms_cfg.ms1_tolerance, ...
                'minMSMSnum', msms_cfg.min_MSMS_num, ...
                'alpha', msms_cfg.alpha, ...
                'resFilterThres', msms_cfg.result_filter_threshold));
            executor = CIMPProcessingExecutor(pipeline_cfg);

            CLogger.info('Re-quantifying at peptide level...');
            rec_counter = 0;
            for idx_block = 1:numel(report.blocks)
                prot_name = '';
                if ~isempty(report.blocks(idx_block).protein_name_pos)
                    prot_name = report.blocks(idx_block).protein_name_pos{1,1};
                end
                records = report.blocks(idx_block).records;
                for idx_rec = 1:numel(records)
                    rec_counter = rec_counter + 1;
                    print_progress = print_progress.update_show(rec_counter);

                    mgf_name = records(idx_rec).raw_name;
                    current_charge = records(idx_rec).charge;
                    current_peptide = records(idx_rec).imp_name;
                    lfMass = get_mass_peptide(current_peptide);
                    lfMz = (lfMass + current_charge * CConstant.pmass) / current_charge;
                    current_key = CIMPQuantRecord.build_id(current_peptide, current_charge, mgf_name);
                    if pep_prot_map.isKey(current_key)
                        prot_name = pep_prot_map(current_key);
                    end
                    rawIdentManager = CIMPRawIdentManager();

                    rt_peaks = records(idx_rec).rt_peaks;
                    for idx_peak = 1:numel(rt_peaks)
                        rt_median = (rt_peaks(idx_peak).rt_start + rt_peaks(idx_peak).rt_end) / 2;
                        rawStore = rawIdentManager.getOrCreate(mgf_name);
                        rawStore = rawStore.appendSpecQuant(rt_median, 1, lfMz, current_charge, {current_peptide}, lfMass, 1);
                        rawIdentManager.setStore(mgf_name, rawStore);
                    end

                    block = executor.requantifyPeptideBlock({prot_name, -1}, rawIdentManager, pep_rtrange_map);
                    report_requant = report_requant.append_block(block);
                end
            end

            print_progress.last_update();
            CLogger.info('Normalization peptide-level re-quantification done.');

            CIMPQuantResultIO.write(report_requant, output_path);
        end
    end
end

function lfMass = get_mass_peptide(pep_seq)
% Input:
%   pep_seq (1 x 1 char/string)
%       peptide sequence
% Output:
%   lfMass (1 x 1 double) Da
%       monoisotopic peptide mass including H2O
lfMass = sum(CConstant.vAAmass(pep_seq-'A'+1));
lfMass = lfMass + CConstant.hmass*2 + CConstant.omass;
end
