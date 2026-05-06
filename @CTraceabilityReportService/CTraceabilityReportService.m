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
            obj.writePeptideMsmsTrace();
            obj.writeSitePeptideTrace();
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
            result = CMS2ResultIO.read(obj.m_cfg.msms_res_path);
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
            for idx_pep = 1:numel(result.Peptides)
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
                    end
                end
            end

            CLogger.info(['[CTraceabilityReportService:writePeptideMsmsTrace] ', ...
                'Trace written: %s'], obj.m_cfg.output_trace_peptide_msms_path);
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

            % The peptide-level report is block-oriented: a protein context
            % line applies to the following IMP records until the next protein
            % context line appears. For traceability, the context is preserved
            % verbatim instead of recalculating protein positions here.
            while ~feof(fin)
                strline = fgetl(fin);
                line_no = line_no + 1;
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
                    end
                elseif strline(1) ~= '@'
                    current_protein_line = strline;
                end
            end

            CLogger.info(['[CTraceabilityReportService:writeSitePeptideTrace] ', ...
                'Trace written: %s'], obj.m_cfg.output_trace_site_peptide_path);
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
                charge_map = obj.m_mgf_dataset_io.get_dataset_charge_map(dataset_name);
                obj.m_dataset_charge_maps(dataset_name) = charge_map;
            end

            if ~isKey(charge_map, spectrum_name)
                error('CTraceabilityReportService:MissingSpectrumCharge', ...
                    ['Cannot find precursor charge for dataset="%s", spectrum_name="%s" ', ...
                    'while reading MSMS result "%s".'], ...
                    dataset_name, spectrum_name, obj.m_cfg.msms_res_path);
            end
            charge = charge_map(spectrum_name);
        end

        function site_records = buildSiteRecords(obj, imp_name, protein_site_positions)
            % Build site contribution records for one IMP sequence.
            % Input:
            %   obj (CTraceabilityReportService)
            %       Service instance with modification abbreviation settings.
            %   imp_name (char/string)
            %       Modified peptide/IMP string from the peptide-level report.
            %   protein_site_positions (char/string)
            %       Protein-site context copied from the peptide-level report.
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

            positions_seq = zeros(1, numel(mod_matches));
            positions_str = zeros(1, numel(mod_matches));
            start_pos = 0;

            % Identify target modification occurrences in the IMP string. The
            % parsed sequence/string positions are used only to recover the
            % modified residue or terminal label; protein positions are not
            % recomputed in this traceability stage.
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
                    else
                        residue_or_term = 'C-term';
                    end
                end
                site_name = [protein_site_positions, ' ', residue_or_term, '_', mod_abbr];

                new_idx = numel(site_records) + 1;
                site_records(new_idx).site_name = site_name;
                site_records(new_idx).protein_site_positions = protein_site_positions;
                site_records(new_idx).residue_or_term = residue_or_term;
                site_records(new_idx).mod_name = mod_name;
                site_records(new_idx).mod_abbr = mod_abbr;
            end
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
                'mod_name_abbr', 'ignore_strings', 'column_idxs'};
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
