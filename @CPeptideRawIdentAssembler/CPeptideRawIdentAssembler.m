classdef CPeptideRawIdentAssembler
    % Build CIMPRawIdentManager instances from peptide-related sources.
    % This class only assembles data and does not execute business workflows.

    methods (Static)
        function [rawIdentManager, filterStats] = buildFromSpectrumList(spectrum_list, deps)
            % Build a raw identification manager from a spectrum list.
            % Input:
            %   spectrum_list (struct array)
            %       spectrum entries with fields:
            %       - dataset_name
            %       - spectrum_name
            %       - peptidoform_list_str
            %       - peptidoform_list_abun
            %   deps (struct)
            %       dependencies:
            %       - getProfilesFunc (function_handle)
            %           function handle: [rt,iso,mz,ch] = f(dataset_name, spectrum_name)
            %           rt: MS1 RT of precursor scan
            %           iso: reference isotope intensity near precursor
            %           mz: precursor m/z
            %           ch: precursor charge
            %       - fixedModNameMass (Nf x 3 cell)
            %           fixed modification list {mod_name, specificity, mass}
            %       - variableModNameMass (Nv x 3 cell)
            %           variable modification list {mod_name, specificity, mass}
            %       - msmsStabilityFilter (optional, 1 x 1 struct)
            %           MSMS stability filter config with fields:
            %           * enabled (logical)
            %           * min_jaccard (1 x 1 double or [])
            %           * min_support_frequency (1 x 1 double or [])
            %           * max_abundance_mad (1 x 1 double or [])
            %           * nan_as_fail (logical)
            % Output:
            %   rawIdentManager (CIMPRawIdentManager)
            %       per-raw identification store manager
            %   filterStats (1 x 1 struct)
            %       filtering counters:
            %       - total_spectra
            %       - kept_spectra
            %       - dropped_spectra_jaccard
            %       - total_imp_candidates
            %       - kept_imp_candidates
            %       - dropped_imp_candidates

            CPeptideRawIdentAssembler.assertDepsForSpectrumList(deps);

            rawIdentManager = CIMPRawIdentManager();
            modNameMass = [deps.fixedModNameMass; deps.variableModNameMass];
            filterCfg = CPeptideRawIdentAssembler.resolveMsmsStabilityFilterCfg(deps);
            filterStats = struct( ...
                'total_spectra', length(spectrum_list), ...
                'kept_spectra', 0, ...
                'dropped_spectra_jaccard', 0, ...
                'total_imp_candidates', 0, ...
                'kept_imp_candidates', 0, ...
                'dropped_imp_candidates', 0);

            for idx_spec = 1:length(spectrum_list)
                spectrum_entry = spectrum_list(idx_spec);
                [is_kept, spectrum_entry, spectrum_stats] = ...
                    CPeptideRawIdentAssembler.applyMsmsStabilityFilterToSpectrum(spectrum_entry, filterCfg);
                filterStats.dropped_spectra_jaccard = filterStats.dropped_spectra_jaccard + spectrum_stats.dropped_spectra_jaccard;
                filterStats.total_imp_candidates = filterStats.total_imp_candidates + spectrum_stats.total_imp_candidates;
                filterStats.kept_imp_candidates = filterStats.kept_imp_candidates + spectrum_stats.kept_imp_candidates;
                filterStats.dropped_imp_candidates = filterStats.dropped_imp_candidates + spectrum_stats.dropped_imp_candidates;

                if ~is_kept
                    continue;
                end

                filterStats.kept_spectra = filterStats.kept_spectra + 1;

                dataset_name = spectrum_entry.dataset_name;
                spectrum_name = spectrum_entry.spectrum_name;
                peptidoform_strs = spectrum_entry.peptidoform_list_str;
                peptidoform_abuns = spectrum_entry.peptidoform_list_abun;

                [isorts, c_ref_isointens, c_mz, cur_ch] = deps.getProfilesFunc(dataset_name, spectrum_name);
                lfMasses = CPeptideRawIdentAssembler.getMassesIMPs(peptidoform_strs, modNameMass);

                rawStore = rawIdentManager.getOrCreate(dataset_name);
                rawStore = rawStore.appendSpecQuant(isorts, c_ref_isointens, c_mz, cur_ch, ...
                    peptidoform_strs, lfMasses, peptidoform_abuns);
                rawIdentManager.setStore(dataset_name, rawStore);
            end
        end

        function pep_quant = buildFromFdrEntries(entries, peptide_list, pep_quant, deps)
            % Parse FDR filtered entries and update per-peptide raw managers.
            % Input:
            %   entries (struct array)
            %       FDR filtered entries
            %   peptide_list (K x 1 cellstr/string)
            %       target peptide sequences
            %   pep_quant (K x 1 cell)
            %       per-peptide raw identification managers
            %   deps (struct)
            %       dependencies:
            %       - ms12DatasetIO (CMS12DatasetIO)
            %           reader providing MS1/MS2 indices and peaks
            %       - ms1_tolerance (struct)
            %           tolerance struct with fields:
            %           * value (1 x 1 double)
            %           * isppm (logical/double)
            % Output:
            %   pep_quant (K x 1 cell)
            %       updated per-peptide managers

            CPeptideRawIdentAssembler.assertDepsForFdrEntries(deps);

            for i_entry = 1:length(entries)
                entry = entries(i_entry);

                if ~isequal(entry.modification, '-')
                    continue;
                end

                for i_list = 1:length(peptide_list)
                    if ~isequal(entry.peptide, peptide_list{i_list})
                        continue;
                    end

                    curr_change = str2double(entry.Charge);
                    curr_mz = str2double(entry.precursor_neutral_mass) / curr_change + CConstant.pmass;
                    curr_MS2_scan = str2double(entry.Scan);

                    mgf_name = erase(entry.DatasetName, '.mgf');
                    ms2_name = deps.ms12DatasetIO.m_cMsFileMapper.get_ms2_stem(mgf_name);
                    MS2_index = deps.ms12DatasetIO.m_mapNameMS2Index(ms2_name);

                    tmp_idx = MS2_index(:,2) == curr_MS2_scan;
                    if ~any(tmp_idx)
                        error('CPeptideRawIdentAssembler:MissingMS2Scan', ...
                            'No matching MS2 scan found for scan number %d in %s', curr_MS2_scan, mgf_name);
                    end

                    MS1Scan = MS2_index(tmp_idx, 1);
                    MS1_index = deps.ms12DatasetIO.m_mapNameMS1Index(ms2_name);
                    MS1_peaks = deps.ms12DatasetIO.m_mapNameMS1Peaks(ms2_name);

                    ino = find(MS1_index(:,1) == MS1Scan);
                    if isempty(ino)
                        error('CPeptideRawIdentAssembler:MissingMS1Scan', ...
                            'No matching MS1 scan found for scan number %d in %s', MS1Scan, mgf_name);
                    end

                    curr_rt = MS1_index(ino, 2);

                    first_peak_idx = [1; MS1_index(1:size(MS1_index), 3)];
                    IX = first_peak_idx(ino):first_peak_idx(ino+1)-1;
                    mz = MS1_peaks(IX, 1);
                    inten = MS1_peaks(IX, 2);

                    if deps.ms1_tolerance.isppm
                        mz_tol = deps.ms1_tolerance.value * curr_mz * 1e-6;
                    else
                        mz_tol = deps.ms1_tolerance.value;
                    end

                    cur_inten = max(inten(mz >= curr_mz-mz_tol & mz <= curr_mz+mz_tol));
                    if isempty(cur_inten)
                        break;
                    end

                    lfMass = CPeptideRawIdentAssembler.getMassPeptide(peptide_list{i_list});
                    rawStore = pep_quant{i_list}.getOrCreate(entry.DatasetName);
                    rawStore = rawStore.appendSpecQuant(curr_rt, cur_inten, curr_mz, curr_change, ...
                        {entry.peptide}, lfMass, 1);
                    pep_quant{i_list}.setStore(entry.DatasetName, rawStore);
                    break;
                end
            end
        end
    end

    methods (Static, Access = private)
        function assertDepsForSpectrumList(deps)
            % Validate dependencies for buildFromSpectrumList.
            % Input:
            %   deps (struct)
            %       dependency struct containing required fields
            required_fields = {'getProfilesFunc', 'fixedModNameMass', 'variableModNameMass'};
            CPeptideRawIdentAssembler.assertRequiredFields(deps, required_fields, 'buildFromSpectrumList');
            if ~isa(deps.getProfilesFunc, 'function_handle')
                error('CPeptideRawIdentAssembler:InvalidDependency', ...
                    'deps.getProfilesFunc must be a function handle.');
            end
        end

        function assertDepsForFdrEntries(deps)
            % Validate dependencies for buildFromFdrEntries.
            % Input:
            %   deps (struct)
            %       dependency struct containing required fields
            required_fields = {'ms12DatasetIO', 'ms1_tolerance'};
            CPeptideRawIdentAssembler.assertRequiredFields(deps, required_fields, 'buildFromFdrEntries');
            if ~isstruct(deps.ms1_tolerance) || ~isfield(deps.ms1_tolerance, 'value') || ~isfield(deps.ms1_tolerance, 'isppm')
                error('CPeptideRawIdentAssembler:InvalidDependency', ...
                    'deps.ms1_tolerance must be a struct with fields: value, isppm.');
            end
        end

        function assertRequiredFields(s, required_fields, method_name)
            % Validate that all required fields exist in a struct.
            % Input:
            %   s (struct)
            %       input struct to validate
            %   required_fields (1 x N cell)
            %       required field names
            %   method_name (1 x 1 char/string)
            %       caller method name for error message
            if nargin < 1 || ~isstruct(s)
                error('CPeptideRawIdentAssembler:InvalidDependency', ...
                    'deps must be a struct for %s.', method_name);
            end
            for idx = 1:numel(required_fields)
                if ~isfield(s, required_fields{idx})
                    error('CPeptideRawIdentAssembler:MissingDependency', ...
                        'Missing deps.%s for %s.', required_fields{idx}, method_name);
                end
            end
        end

        function lfMasses = getMassesIMPs(cstrIMP, modNameMass)
            % Get masses of each IMP.
            % Input:
            %   cstrIMP (1 x N cell)
            %       IMP sequence strings with modification tags
            %   modNameMass (M x 3 cell)
            %       modification list {mod_name, specificity, mass}
            % Output:
            %   lfMasses (N x 1 double)
            %       calculated monoisotopic masses including H2O
            lfMasses = zeros(length(cstrIMP),1);
            for idx_imp = 1:length(cstrIMP)
                mod_seq = cstrIMP{idx_imp};
                reg_exp = '\{(.*?)\}';
                [mod_str, seq_str] = regexp(mod_seq, reg_exp, 'tokens', 'split');

                seq_str = strjoin(seq_str, '');
                seq_str([1,end]) = [];
                lfMasses(idx_imp) = sum(CConstant.vAAmass(seq_str-'A'+1));

                for idx_mod = 1:length(mod_str)
                    is_notfound = true;
                    for idx_mlist = 1:size(modNameMass,1)
                        if isequal(modNameMass{idx_mlist,1}, mod_str{idx_mod}{1})
                            is_notfound = false;
                            lfMasses(idx_imp) = lfMasses(idx_imp) + modNameMass{idx_mlist,3};
                            break;
                        end
                    end
                    if is_notfound
                        error('CPeptideRawIdentAssembler:UnexpectedModification', ...
                            'Unexpected modification is found: "%s".', mod_str{idx_mod}{1});
                    end
                end

                lfMasses(idx_imp) = lfMasses(idx_imp) + CConstant.hmass*2 + CConstant.omass;
            end
        end

        function lfMass = getMassPeptide(pep_seq)
            % Get mass of peptide sequence.
            % Input:
            %   pep_seq (1 x 1 char/string)
            %       peptide sequence
            % Output:
            %   lfMass (1 x 1 double)
            %       monoisotopic peptide mass including H2O
            lfMass = sum(CConstant.vAAmass(pep_seq-'A'+1));
            lfMass = lfMass + CConstant.hmass*2 + CConstant.omass;
        end

        function filterCfg = resolveMsmsStabilityFilterCfg(deps)
            % Resolve optional MSMS stability filter config.
            % Input:
            %   deps (struct)
            %       dependency struct potentially containing msmsStabilityFilter
            % Output:
            %   filterCfg (1 x 1 struct)
            %       normalized filter config
            filterCfg = struct( ...
                'enabled', false, ...
                'min_jaccard', [], ...
                'min_support_frequency', [], ...
                'max_abundance_mad', [], ...
                'nan_as_fail', true);
            if isfield(deps, 'msmsStabilityFilter') && isstruct(deps.msmsStabilityFilter)
                cfg = deps.msmsStabilityFilter;
                if isfield(cfg, 'enabled')
                    filterCfg.enabled = logical(cfg.enabled);
                end
                if isfield(cfg, 'min_jaccard')
                    filterCfg.min_jaccard = cfg.min_jaccard;
                end
                if isfield(cfg, 'min_support_frequency')
                    filterCfg.min_support_frequency = cfg.min_support_frequency;
                end
                if isfield(cfg, 'max_abundance_mad')
                    filterCfg.max_abundance_mad = cfg.max_abundance_mad;
                end
                if isfield(cfg, 'nan_as_fail')
                    filterCfg.nan_as_fail = logical(cfg.nan_as_fail);
                end
            end
        end

        function [is_kept, filtered_spectrum, stats] = applyMsmsStabilityFilterToSpectrum(spectrum_entry, filterCfg)
            % Apply spectrum-level and IMP-level stability filtering.
            % Input:
            %   spectrum_entry (1 x 1 struct)
            %       one spectrum entry from CMS2Result.Peptides(i).spectrum_list
            %   filterCfg (1 x 1 struct)
            %       normalized filter config
            % Output:
            %   is_kept (logical)
            %       whether this spectrum should be kept
            %   filtered_spectrum (1 x 1 struct)
            %       filtered spectrum entry
            %   stats (1 x 1 struct)
            %       filtering counters for this spectrum
            filtered_spectrum = spectrum_entry;
            total_imp = numel(spectrum_entry.peptidoform_list_str);
            stats = struct( ...
                'dropped_spectra_jaccard', 0, ...
                'total_imp_candidates', total_imp, ...
                'kept_imp_candidates', total_imp, ...
                'dropped_imp_candidates', 0);

            if ~filterCfg.enabled
                is_kept = true;
                return;
            end

            if ~CPeptideRawIdentAssembler.passesMinMetric(spectrum_entry, 'jaccard_stability', ...
                    filterCfg.min_jaccard, filterCfg.nan_as_fail)
                is_kept = false;
                stats.dropped_spectra_jaccard = 1;
                stats.kept_imp_candidates = 0;
                stats.dropped_imp_candidates = total_imp;
                return;
            end

            keep_mask = true(1, total_imp);
            support_vals = CPeptideRawIdentAssembler.getOptionalNumericArray(spectrum_entry, ...
                'peptidoform_list_support_freq', total_imp);
            mad_vals = CPeptideRawIdentAssembler.getOptionalNumericArray(spectrum_entry, ...
                'peptidoform_list_abundance_mad', total_imp);

            if ~isempty(filterCfg.min_support_frequency)
                keep_mask = keep_mask & CPeptideRawIdentAssembler.passesMinArray( ...
                    support_vals, filterCfg.min_support_frequency, filterCfg.nan_as_fail);
            end
            if ~isempty(filterCfg.max_abundance_mad)
                keep_mask = keep_mask & CPeptideRawIdentAssembler.passesMaxArray( ...
                    mad_vals, filterCfg.max_abundance_mad, filterCfg.nan_as_fail);
            end

            filtered_spectrum.peptidoform_list_str = spectrum_entry.peptidoform_list_str(keep_mask);
            filtered_spectrum.peptidoform_list_abun = spectrum_entry.peptidoform_list_abun(keep_mask);
            if isfield(spectrum_entry, 'peptidoform_list_support_freq')
                filtered_spectrum.peptidoform_list_support_freq = support_vals(keep_mask);
            end
            if isfield(spectrum_entry, 'peptidoform_list_abundance_mad')
                filtered_spectrum.peptidoform_list_abundance_mad = mad_vals(keep_mask);
            end

            stats.kept_imp_candidates = sum(keep_mask);
            stats.dropped_imp_candidates = total_imp - stats.kept_imp_candidates;
            is_kept = ~isempty(filtered_spectrum.peptidoform_list_str);
        end

        function pass_mask = passesMinArray(values, threshold, nan_as_fail)
            % Evaluate array-wise lower-bound rule.
            pass_mask = values >= threshold;
            if nan_as_fail
                pass_mask(isnan(values)) = false;
            else
                pass_mask(isnan(values)) = true;
            end
        end

        function pass_mask = passesMaxArray(values, threshold, nan_as_fail)
            % Evaluate array-wise upper-bound rule.
            pass_mask = values <= threshold;
            if nan_as_fail
                pass_mask(isnan(values)) = false;
            else
                pass_mask(isnan(values)) = true;
            end
        end

        function pass = passesMinMetric(spectrum_entry, field_name, threshold, nan_as_fail)
            % Evaluate scalar lower-bound rule from optional spectrum field.
            if isempty(threshold)
                pass = true;
                return;
            end

            if ~isfield(spectrum_entry, field_name) || isempty(spectrum_entry.(field_name))
                pass = ~nan_as_fail;
                return;
            end

            metric_val = spectrum_entry.(field_name);
            if isnan(metric_val)
                pass = ~nan_as_fail;
                return;
            end
            pass = metric_val >= threshold;
        end

        function values = getOptionalNumericArray(spectrum_entry, field_name, expected_len)
            % Get optional numeric array with safe default NaN padding.
            values = NaN(1, expected_len);
            if ~isfield(spectrum_entry, field_name)
                return;
            end

            raw_values = spectrum_entry.(field_name);
            if isempty(raw_values)
                return;
            end

            take_len = min(numel(raw_values), expected_len);
            values(1:take_len) = raw_values(1:take_len);
        end
    end
end
