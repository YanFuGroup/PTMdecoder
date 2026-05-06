classdef CTraceabilityReportService < handle
    % Generate traceability sidecar reports without changing primary outputs.
    %
    % This service writes two tab-separated files:
    %   1. report_trace_peptide_msms.txt links charged IMP entries to their
    %      supporting MS/MS spectra.
    %   2. report_trace_site_peptide.txt links site-level dataset-matrix keys
    %      to peptide-level IMP quantification records.
    %
    % The service is intentionally implemented as a sidecar writer. It reads
    % existing MS/MS and peptide-level reports, preserves their public formats,
    % and emits additional long-table mappings for evidence tracing.

    properties (Access = private)
        m_cfg
        m_mgf_dataset_io
        m_dataset_charge_maps
    end

    methods
        function obj = CTraceabilityReportService(cfg)
            % Construct traceability report service.
            % Input:
            %   cfg (struct)
            %       Normalized config from CTraceabilityReportServiceConfig.
            if ~isstruct(cfg)
                error('CTraceabilityReportService:InvalidConstructorArgs', ...
                    'Expected a config struct.');
            end
            obj.m_cfg = obj.finalizeConfig(cfg);
            obj.m_mgf_dataset_io = CMgfDatasetIO(obj.m_cfg.spec_dir_path);
            obj.m_dataset_charge_maps = containers.Map('KeyType', 'char', 'ValueType', 'any');
        end

        function run(obj)
            % Generate all traceability sidecar reports.
            % Input:
            %   obj (CTraceabilityReportService)
            %       Configured service instance.
            % Output:
            %   (none)
            %       Writes the configured peptide-MS/MS and site-peptide trace
            %       files to disk.
            CLogger.info('[CTraceabilityReportService:run] Traceability report generation started.');
            total_timer = tic;

            CLogger.info('[CTraceabilityReportService:run] Stage start: peptide-msms trace.');
            stage_timer = tic;
            obj.writePeptideMsmsTrace();
            CLogger.info('[CTraceabilityReportService:run] Stage end: peptide-msms trace. elapsed=%.1fs', ...
                toc(stage_timer));

            CLogger.info('[CTraceabilityReportService:run] Stage start: site-peptide trace.');
            stage_timer = tic;
            obj.writeSitePeptideTrace();
            CLogger.info('[CTraceabilityReportService:run] Stage end: site-peptide trace. elapsed=%.1fs', ...
                toc(stage_timer));

            CLogger.info('[CTraceabilityReportService:run] Traceability report generation done. elapsed=%.1fs', ...
                toc(total_timer));
        end

        function delete(obj)
            % Release opened MGF file handles owned by this service.
            if ~isempty(obj.m_mgf_dataset_io)
                delete(obj.m_mgf_dataset_io);
            end
        end
    end

    methods (Access = private)
        function writePeptideMsmsTrace(obj)
            % Write charged IMP-to-MSMS evidence trace.
            % Input:
            %   obj (CTraceabilityReportService)
            %       Configured service instance with MS/MS report path, MGF
            %       directory, and output path.
            % Output:
            %   (none)
            %       Writes one row per peptide/spectrum/IMP evidence relation.
            CLogger.info('[CTraceabilityReportService:writePeptideMsmsTrace] Reading MS/MS result: %s', ...
                obj.m_cfg.msms_res_path);
            result = CMS2ResultIO.read(obj.m_cfg.msms_res_path);
            CLogger.info('[CTraceabilityReportService:writePeptideMsmsTrace] MS/MS result loaded. peptides=%d', ...
                numel(result.Peptides));
            obj.ensureParentDirForFile(obj.m_cfg.output_trace_peptide_msms_path);

            fid = fopen(obj.m_cfg.output_trace_peptide_msms_path, 'w');
            if fid < 0
                error('CTraceabilityReportService:OpenPeptideMsmsTraceFailed', ...
                    'Cannot open output file: %s', obj.m_cfg.output_trace_peptide_msms_path);
            end
            cleanup = onCleanup(@() fclose(fid));

            fprintf(fid, ['imp_msms_key\tmsms_key\tdataset_name\timp_name\tcharge\t', ...
                'peptide_sequence\tspectrum_name\tmsms_abundance\tsupport_frequency\tabundance_mad\n']);

            % The MS/MS report does not store precursor charge. Charge is
            % therefore resolved from the corresponding MGF dataset by exact
            % spectrum name. A missing charge is treated as an input mismatch
            % rather than silently falling back to a scan-number lookup.
            print_progress = CPrintProgress(max(numel(result.Peptides), 1), 'trace_peptide_msms');
            row_count = 0;
            for idx_pep = 1:numel(result.Peptides)
                print_progress = print_progress.update_show(idx_pep);
                peptide = result.Peptides(idx_pep);
                for idx_spec = 1:numel(peptide.spectrum_list)
                    spectrum = peptide.spectrum_list(idx_spec);
                    dataset_name = char(string(spectrum.dataset_name));
                    spectrum_name = char(string(spectrum.spectrum_name));
                    charge = obj.lookupSpectrumCharge(dataset_name, spectrum_name);
                    charge_str = num2str(charge);
                    msms_key = obj.joinKey({dataset_name, spectrum_name});

                    for idx_imp = 1:spectrum.peptidoform_num
                        imp_name = char(string(spectrum.peptidoform_list_str{idx_imp}));
                        imp_msms_key = obj.joinKey({dataset_name, imp_name, ['+', charge_str]});
                        abundance = spectrum.peptidoform_list_abun(idx_imp);
                        support_frequency = obj.getVectorValueOrNaN( ...
                            spectrum, 'peptidoform_list_support_freq', idx_imp);
                        abundance_mad = obj.getVectorValueOrNaN( ...
                            spectrum, 'peptidoform_list_abundance_mad', idx_imp);

                        fprintf(fid, '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%.15g\t%.15g\t%.15g\n', ...
                            imp_msms_key, ...
                            msms_key, ...
                            dataset_name, ...
                            imp_name, ...
                            charge_str, ...
                            char(string(peptide.peptide_sequence)), ...
                            spectrum_name, ...
                            abundance, ...
                            support_frequency, ...
                            abundance_mad);
                        row_count = row_count + 1;
                    end
                end
            end
            print_progress.last_update();

            CLogger.info(['[CTraceabilityReportService:writePeptideMsmsTrace] ', ...
                'Trace written: %s, rows=%d'], obj.m_cfg.output_trace_peptide_msms_path, row_count);
        end

        function writeSitePeptideTrace(obj)
            % Write site-to-peptide/IMP contribution trace.
            % Input:
            %   obj (CTraceabilityReportService)
            %       Configured service instance with peptide-level report path,
            %       site-modification mapping, and output path.
            % Output:
            %   (none)
            %       Writes one row per site/IMP quantification contribution.
            obj.ensureParentDirForFile(obj.m_cfg.output_trace_site_peptide_path);

            fin = fopen(obj.m_cfg.pep_level_file_path, 'r');
            if fin < 0
                error('CTraceabilityReportService:OpenPeptideInputFailed', ...
                    'Cannot open peptide-level result file: %s', obj.m_cfg.pep_level_file_path);
            end
            cleanup_in = onCleanup(@() fclose(fin));
            pep_file_info = dir(obj.m_cfg.pep_level_file_path);
            pep_file_size = 1;
            if ~isempty(pep_file_info)
                pep_file_size = max(pep_file_info.bytes, 1);
            end

            fout = fopen(obj.m_cfg.output_trace_site_peptide_path, 'w');
            if fout < 0
                error('CTraceabilityReportService:OpenSitePeptideTraceFailed', ...
                    'Cannot open output file: %s', obj.m_cfg.output_trace_site_peptide_path);
            end
            cleanup_out = onCleanup(@() fclose(fout));

            fprintf(fout, ['site_peptide_key\tsite_name\tdataset_name\timp_name\tcharge\t', ...
                'peptide_sequence\tprotein_context\tprotein_site_positions\tresidue_or_term\t', ...
                'mod_name\tmod_abbr\tmass_center\tlow_mz_bound\thigh_mz_bound\tpeak_area\tpeptide_line_no\n']);

            % Keep compatibility with peptide-level report format:
            % line 1: protein header, line 2: peptide header, line 3: xic peak header
            for idx_header = 1:3
                if ~feof(fin)
                    fgetl(fin);
                end
            end

            line_no = 3;
            current_protein_line = '';
            row_count = 0;

            % The peptide-level report is block-oriented: a protein context
            % line applies to the following IMP records until the next protein
            % context line appears. The original context is preserved verbatim,
            % while protein_site_positions stores computed modification sites.
            while ~feof(fin)
                strline = fgetl(fin);
                line_no = line_no + 1;
                if mod(line_no, 2000) == 0
                    CLogger.progress('trace_site_peptide', ftell(fin), pep_file_size);
                end
                if ~ischar(strline) || isempty(strline)
                    continue;
                end

                if strline(1) == '*'
                    rec = obj.parsePeptideRecordLine(strline);
                    site_records = obj.buildSiteRecords(rec.imp_name, current_protein_line);
                    peptide_sequence = obj.buildPlainPeptideSequence(rec.imp_name);
                    for idx_site = 1:numel(site_records)
                        site_rec = site_records(idx_site);
                        charge_str = num2str(rec.charge);
                        site_peptide_key = obj.joinKey({site_rec.site_name, rec.dataset_name, ...
                            rec.imp_name, ['+', charge_str]});

                        fprintf(fout, '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%.15g\t%.15g\t%.15g\t%.15g\t%d\n', ...
                            site_peptide_key, ...
                            site_rec.site_name, ...
                            rec.dataset_name, ...
                            rec.imp_name, ...
                            charge_str, ...
                            peptide_sequence, ...
                            current_protein_line, ...
                            site_rec.protein_site_positions, ...
                            site_rec.residue_or_term, ...
                            site_rec.mod_name, ...
                            site_rec.mod_abbr, ...
                            rec.mass_center, ...
                            rec.low_mz_bound, ...
                            rec.high_mz_bound, ...
                            rec.peak_area, ...
                            line_no);
                        row_count = row_count + 1;
                    end
                elseif strline(1) ~= '@'
                    current_protein_line = strline;
                end
            end
            CLogger.progress('trace_site_peptide', pep_file_size, pep_file_size);

            CLogger.info(['[CTraceabilityReportService:writeSitePeptideTrace] ', ...
                'Trace written: %s, rows=%d'], obj.m_cfg.output_trace_site_peptide_path, row_count);
        end

        function charge = lookupSpectrumCharge(obj, dataset_name, spectrum_name)
            % Look up spectrum precursor charge by exact spectrum name.
            % Input:
            %   obj (CTraceabilityReportService)
            %       Service instance that owns the MGF dataset reader.
            %   dataset_name (char/string)
            %       MGF dataset filename from report_msms.txt.
            %   spectrum_name (char/string)
            %       Spectrum name from report_msms.txt.
            % Output:
            %   charge (double)
            %       Precursor charge read from the MGF TITLE/CHARGE map.
            if isKey(obj.m_dataset_charge_maps, dataset_name)
                charge_map = obj.m_dataset_charge_maps(dataset_name);
            else
                CLogger.info('[CTraceabilityReportService:lookupSpectrumCharge] Loading charge map: %s', ...
                    dataset_name);
                charge_map = obj.m_mgf_dataset_io.get_dataset_charge_map(dataset_name);
                obj.m_dataset_charge_maps(dataset_name) = charge_map;
                CLogger.info('[CTraceabilityReportService:lookupSpectrumCharge] Charge map ready: %s, entries=%d', ...
                    dataset_name, charge_map.Count);
            end

            if ~isKey(charge_map, spectrum_name)
                error('CTraceabilityReportService:MissingSpectrumCharge', ...
                    ['Cannot find precursor charge for dataset="%s", spectrum_name="%s" ', ...
                    'while reading MSMS result "%s".'], ...
                    dataset_name, spectrum_name, obj.m_cfg.msms_res_path);
            end
            charge = charge_map(spectrum_name);
        end

        function site_records = buildSiteRecords(obj, imp_name, protein_context)
            % Build site contribution records for one IMP sequence.
            % Input:
            %   obj (CTraceabilityReportService)
            %       Service instance with modification abbreviation settings.
            %   imp_name (char/string)
            %       Modified peptide/IMP string from the peptide-level report.
            %   protein_context (char/string)
            %       Protein/start-position context copied from the peptide-level report.
            % Output:
            %   site_records (struct array)
            %       One record per target modification occurrence in imp_name.
            modified_peptide = char(string(imp_name));
            for idx_ig = 1:numel(obj.m_cfg.ignore_strings)
                modified_peptide = strrep(modified_peptide, obj.m_cfg.ignore_strings{idx_ig}, '');
            end

            mod_matches = regexp(modified_peptide, '{(.*?)}', 'tokens');
            site_records = struct('site_name', {}, 'protein_site_positions', {}, 'residue_or_term', {}, ...
                'mod_name', {}, 'mod_abbr', {});
            if isempty(mod_matches)
                return;
            end

            protein_contexts = obj.parseProteinContexts(protein_context);
            if isempty(protein_contexts)
                return;
            end
            peptide_length = numel(obj.buildPlainPeptideSequence(modified_peptide));

            positions_seq = zeros(1, numel(mod_matches));
            positions_str = zeros(1, numel(mod_matches));
            start_pos = 0;

            % Identify target modification occurrences in the IMP string and
            % map each occurrence back to every protein context for tracing.
            for idx_mod = 1:numel(mod_matches)
                mod_name = mod_matches{idx_mod}{1};
                found_index = strfind(modified_peptide(start_pos + 1:end), ['{' mod_name '}']);
                if isempty(found_index)
                    continue;
                end

                if idx_mod == 1
                    positions_seq(idx_mod) = found_index(1) - 2;
                else
                    positions_seq(idx_mod) = positions_seq(idx_mod - 1) + found_index(1) - 1;
                end

                positions_str(idx_mod) = start_pos + found_index(1);
                start_pos = positions_str(idx_mod) + numel(mod_name) + 1;

                if ~isKey(obj.m_cfg.mod_name_abbr, mod_name)
                    continue;
                end

                mod_abbr = obj.m_cfg.mod_name_abbr(mod_name);
                residue_or_term = modified_peptide(positions_str(idx_mod) - 1);

                if residue_or_term == '_'
                    if positions_seq(idx_mod) == 0
                        residue_or_term = 'N-term';
                        local_pos = 0;
                    else
                        residue_or_term = 'C-term';
                        local_pos = peptide_length + 1;
                    end
                else
                    local_pos = positions_seq(idx_mod);
                end
                protein_site_positions = obj.buildProteinSitePositions(protein_contexts, local_pos);
                site_name = [protein_site_positions, ' ', residue_or_term, '_', mod_abbr];

                new_idx = numel(site_records) + 1;
                site_records(new_idx).site_name = site_name;
                site_records(new_idx).protein_site_positions = protein_site_positions;
                site_records(new_idx).residue_or_term = residue_or_term;
                site_records(new_idx).mod_name = mod_name;
                site_records(new_idx).mod_abbr = mod_abbr;
            end
        end

        function protein_contexts = parseProteinContexts(~, protein_context)
            % Parse protein/start-position pairs from a peptide-level context line.
            % Input:
            %   protein_context (char/string)
            %       Context in the form protein,start;protein,start;.
            % Output:
            %   protein_contexts (struct array)
            %       Valid protein/start-position entries from the context.
            protein_contexts = struct('protein_name', {}, 'start_pos', {});
            segments = strsplit(char(string(protein_context)), ';');
            for idx_seg = 1:(numel(segments) - 1)
                key_value = strsplit(segments{idx_seg}, ',');
                if numel(key_value) < 2
                    continue;
                end
                protein_name = strtrim(key_value{1});
                start_pos = str2double(strtrim(key_value{2}));
                if isempty(protein_name) || isnan(start_pos) || ~isfinite(start_pos)
                    continue;
                end

                new_idx = numel(protein_contexts) + 1;
                protein_contexts(new_idx).protein_name = protein_name;
                protein_contexts(new_idx).start_pos = start_pos;
            end
        end

        function protein_site_positions = buildProteinSitePositions(obj, protein_contexts, local_pos)
            % Build protein-position prefix for one peptide-local site position.
            % Input:
            %   obj (CTraceabilityReportService)
            %       Service instance with site-position numbering config.
            %   protein_contexts (struct array)
            %       Parsed protein/start-position entries.
            %   local_pos (double)
            %       Peptide-local modification position. N-term is 0 and C-term
            %       is peptide length plus one.
            % Output:
            %   protein_site_positions (char)
            %       Protein-site prefix in the form protein,pos;protein,pos;.
            protein_site_tokens = cell(1, numel(protein_contexts));
            for idx_ctx = 1:numel(protein_contexts)
                site_pos = protein_contexts(idx_ctx).start_pos + local_pos - 1;
                if ~obj.m_cfg.site_position_count_initial_m
                    site_pos = site_pos - 1;
                end
                protein_site_tokens{idx_ctx} = [protein_contexts(idx_ctx).protein_name, ...
                    ',', num2str(site_pos), ';'];
            end
            protein_site_positions = strjoin(protein_site_tokens, '');
        end

        function peptide_sequence = buildPlainPeptideSequence(obj, imp_name)
            % Remove modification annotations and terminal markers from an IMP string.
            % Input:
            %   obj (CTraceabilityReportService)
            %       Service instance with ignore-string settings.
            %   imp_name (char/string)
            %       Modified peptide/IMP string.
            % Output:
            %   peptide_sequence (char)
            %       Plain peptide sequence without modification annotations.
            peptide_sequence = char(string(imp_name));
            for idx_ig = 1:numel(obj.m_cfg.ignore_strings)
                peptide_sequence = strrep(peptide_sequence, obj.m_cfg.ignore_strings{idx_ig}, '');
            end
            peptide_sequence = regexprep(peptide_sequence, '\{.*?\}', '');
            peptide_sequence = strrep(peptide_sequence, '_', '');
        end

        function ensureParentDirForFile(~, file_path)
            % Ensure parent directory exists for a target file path.
            % Input:
            %   file_path (char/string)
            %       Output file path whose parent directory should exist.
            % Output:
            %   (none)
            parent_dir = fileparts(char(string(file_path)));
            if isempty(parent_dir)
                return;
            end
            CPathResolver.ensureDir(parent_dir);
        end

        function cfg = finalizeConfig(~, cfg)
            % Validate required config fields.
            % Input:
            %   cfg (struct)
            %       Traceability service configuration.
            % Output:
            %   cfg (struct)
            %       Validated configuration with normalized ignore_strings.
            required_fields = {'spec_dir_path', 'msms_res_path', 'pep_level_file_path', ...
                'output_trace_peptide_msms_path', 'output_trace_site_peptide_path', ...
                'mod_name_abbr', 'ignore_strings', 'site_position_count_initial_m', 'column_idxs'};
            for idx_field = 1:numel(required_fields)
                if ~isfield(cfg, required_fields{idx_field})
                    error('CTraceabilityReportService:MissingConfigField', ...
                        'Missing config field: %s', required_fields{idx_field});
                end
            end
            if ~isa(cfg.mod_name_abbr, 'containers.Map')
                error('CTraceabilityReportService:InvalidModNameAbbrMap', ...
                    'mod_name_abbr must be a containers.Map.');
            end
            if ~iscell(cfg.ignore_strings)
                cfg.ignore_strings = {cfg.ignore_strings};
            end
            if ~islogical(cfg.site_position_count_initial_m) && ...
                    ~(isnumeric(cfg.site_position_count_initial_m) && isscalar(cfg.site_position_count_initial_m))
                error('CTraceabilityReportService:InvalidSitePositionCountInitialM', ...
                    'site_position_count_initial_m must be logical or numeric scalar.');
            end
            cfg.site_position_count_initial_m = logical(cfg.site_position_count_initial_m);
            required_column_fields = {'icol_seq', 'icol_charge', 'icol_dataset', ...
                'icol_mass_center', 'icol_low_mz_bound', 'icol_high_mz_bound', 'icol_auc'};
            for idx_field = 1:numel(required_column_fields)
                if ~isfield(cfg.column_idxs, required_column_fields{idx_field})
                    error('CTraceabilityReportService:InvalidColumnIdxs', ...
                        'column_idxs must include field: %s', required_column_fields{idx_field});
                end
            end
        end
    end

    methods (Access = private)
        function rec = parsePeptideRecordLine(obj, strline)
            % Parse a peptide-level IMP record line.
            % Input:
            %   obj (CTraceabilityReportService)
            %       Service instance with peptide-level column indexes.
            %   strline (char)
            %       Tab-separated peptide-level record line beginning with '*'.
            % Output:
            %   rec (struct)
            %       Parsed IMP quantification fields used by the trace writer.
            segments = strsplit(strline, sprintf('\t'));
            max_column_idx = max([obj.m_cfg.column_idxs.icol_seq, ...
                obj.m_cfg.column_idxs.icol_charge, ...
                obj.m_cfg.column_idxs.icol_dataset, ...
                obj.m_cfg.column_idxs.icol_mass_center, ...
                obj.m_cfg.column_idxs.icol_low_mz_bound, ...
                obj.m_cfg.column_idxs.icol_high_mz_bound, ...
                obj.m_cfg.column_idxs.icol_auc]);
            if numel(segments) < max_column_idx || ~strcmp(segments{1}, '*')
                error('CTraceabilityReportService:InvalidPeptideRecordLine', ...
                    'Invalid peptide-level record line: %s', strline);
            end

            idxs = obj.m_cfg.column_idxs;
            charge = str2double(strrep(strtrim(segments{idxs.icol_charge}), '+', ''));
            mass_center = str2double(strtrim(segments{idxs.icol_mass_center}));
            low_mz_bound = str2double(strtrim(segments{idxs.icol_low_mz_bound}));
            high_mz_bound = str2double(strtrim(segments{idxs.icol_high_mz_bound}));
            peak_area = str2double(strtrim(segments{idxs.icol_auc}));
            if any(isnan([charge, mass_center, low_mz_bound, high_mz_bound, peak_area]))
                error('CTraceabilityReportService:InvalidPeptideRecordNumber', ...
                    'Invalid numeric field in peptide-level record line: %s', strline);
            end

            rec = struct();
            rec.imp_name = strtrim(segments{idxs.icol_seq});
            rec.charge = charge;
            rec.dataset_name = strtrim(segments{idxs.icol_dataset});
            rec.mass_center = mass_center;
            rec.low_mz_bound = low_mz_bound;
            rec.high_mz_bound = high_mz_bound;
            rec.peak_area = peak_area;
        end
    end

    methods (Static, Access = private)
        function value = getVectorValueOrNaN(source_struct, field_name, idx)
            % Safely read one vector value from a spectrum struct.
            % Input:
            %   source_struct (struct)
            %       Spectrum struct parsed from the MS/MS report.
            %   field_name (char/string)
            %       Vector field name to read.
            %   idx (double)
            %       One-based vector index.
            % Output:
            %   value (double)
            %       Field value, or NaN when the field/index is unavailable.
            value = NaN;
            if ~isfield(source_struct, field_name)
                return;
            end
            values = source_struct.(field_name);
            if numel(values) < idx || isempty(values(idx))
                return;
            end
            value = values(idx);
        end

        function key = joinKey(parts)
            % Join stable trace key parts.
            % Input:
            %   parts (cell)
            %       Key components.
            % Output:
            %   key (char)
            %       Pipe-delimited trace key.
            key = strjoin(cellfun(@(x) char(string(x)), parts, 'UniformOutput', false), '|');
        end
    end
end
