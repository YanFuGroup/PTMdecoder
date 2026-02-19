function pep_quant = readFdrPeptides(obj, input_file_path, ms12DatasetIO, peptide_list, pep_quant)
% Parse filtered PSM results and build raw managers for target peptides
% Input:
%   obj (CMSMSPepDeconv)
%       Processor instance
%       File identifier of the filtered result file
%   input_file_path (1 x 1 char/string)
%       Path to the filtered result file
%   ms12DatasetIO (object)
%       MS1/MS2 dataset IO instance
%   peptide_list (1 x K cell)
%       Target peptide sequences
%   pep_quant (1 x K cell)
%       Per-peptide raw identification store managers
% Output:
%   pep_quant (1 x K cell)
%       Updated raw identification managers

if nargin < 3 || isempty(ms12DatasetIO)
    [obj, ~] = obj.ensureMs12DatasetIO();
    ms12DatasetIO = obj.m_cMs12DatasetIO;
end

FDRfilteredResults = CFdrFilteredResultIO.read(input_file_path);
entries = FDRfilteredResults.entries;

for i_entry = 1:length(entries)
    entry = entries(i_entry);

    % Only process peptides without modifications (12th column is '-')
    if ~isequal(entry.modification, '-')
        continue;
    end

    % Check if this peptide is in our target list
    for i_list = 1:length(peptide_list)
        if isequal(entry.peptide, peptide_list{i_list})
            % Extract information from the segment
            curr_change = str2double(entry.Charge);
            curr_mz = str2double(entry.precursor_neutral_mass) / curr_change + CConstant.pmass;
            curr_MS2_scan = str2double(entry.Scan);

            % Find MS2 index
            mgf_name = erase(entry.DatasetName, '.mgf');
            ms2_name = ms12DatasetIO.m_cMsFileMapper.get_ms2_stem(mgf_name);
            MS2_index = ms12DatasetIO.m_mapNameMS2Index(ms2_name);

            % Find the corresponding index of ms2
            tmp_idx = MS2_index(:,2) == curr_MS2_scan;
            if ~any(tmp_idx)
                error('No matching MS2 scan found for scan number %d in %s', curr_MS2_scan, mgf_name);
            end

            % Get MS1 scan number
            MS1Scan = MS2_index(tmp_idx, 1);

            % Get MS1 data
            MS1_index = ms12DatasetIO.m_mapNameMS1Index(ms2_name);
            MS1_peaks = ms12DatasetIO.m_mapNameMS1Peaks(ms2_name);

            ino = find(MS1_index(:,1) == MS1Scan);
            if isempty(ino)
                error('No matching MS1 scan found for scan number %d in %s', MS1Scan, mgf_name);
            end

            curr_rt = MS1_index(ino, 2);

            % Get peak list for this MS1
            first_peak_idx = [1; MS1_index(1:size(MS1_index), 3)];
            IX = first_peak_idx(ino):first_peak_idx(ino+1)-1;
            mz = MS1_peaks(IX, 1);
            inten = MS1_peaks(IX, 2);

            % Calculate tolerance
            if obj.m_ms1_tolerance.isppm
                mz_tol = obj.m_ms1_tolerance.value * curr_mz * 1e-6;
            else
                mz_tol = obj.m_ms1_tolerance.value;
            end

            % Find intensity within tolerance
            cur_inten = max(inten(mz >= curr_mz-mz_tol & mz <= curr_mz+mz_tol));
            if isempty(cur_inten)
                break;
            end

            % Calculate peptide mass
            lfMass = get_mass_peptide(peptide_list{i_list});

            % Add to quantification
            rawStore = pep_quant{i_list}.getOrCreate(entry.DatasetName);
            rawStore = rawStore.appendSpecQuant(curr_rt, cur_inten, curr_mz, curr_change, {entry.peptide}, lfMass, 1);
            pep_quant{i_list}.setStore(entry.DatasetName, rawStore);
            break;
        end
    end
end
end



function lfMass = get_mass_peptide(pep_seq)
% Get the mass of each peptide
lfMass = sum(CConstant.vAAmass(pep_seq-'A'+1));
lfMass = lfMass + CConstant.hmass*2 + CConstant.omass;
end
