classdef CPeptideRawIdentAssembler
    % Build CIMPRawIdentManager instances from peptide-related sources.
    % This class only assembles data and does not execute business workflows.

    methods (Static)
        function rawIdentManager = buildFromSpectrumList(spectrum_list, deps)
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
            % Output:
            %   rawIdentManager (CIMPRawIdentManager)
            %       per-raw identification store manager

            CPeptideRawIdentAssembler.assertDepsForSpectrumList(deps);

            rawIdentManager = CIMPRawIdentManager();
            modNameMass = [deps.fixedModNameMass; deps.variableModNameMass];

            for idx_spec = 1:length(spectrum_list)
                dataset_name = spectrum_list(idx_spec).dataset_name;
                spectrum_name = spectrum_list(idx_spec).spectrum_name;
                peptidoform_strs = spectrum_list(idx_spec).peptidoform_list_str;
                peptidoform_abuns = spectrum_list(idx_spec).peptidoform_list_abun;

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
    end
end
