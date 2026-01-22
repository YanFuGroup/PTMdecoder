function pep_quant = readSearchResult(obj, fin, input_file_path, ms12DatasetIO, pep_quant)
% Process the filtered results file and extract peptide information
    
    progress_printer = CPrintProgress(dir(input_file_path).bytes);
    fgetl(fin); % skip the first line
    
    while ~feof(fin)
        strline = fgetl(fin);
        progress_printer = progress_printer.update_show(ftell(fin));
        
        if isempty(strline)
            continue;
        end

        % Parse the line
        segment = split(strline);
        
        % Only process peptides without modifications (12th column is '-')
        if ~isequal(segment{12}, '-')
            continue;
        end
        
        % Check if this peptide is in our target list
        for i_list = 1:length(obj.peptide_list)
            if isequal(segment{10}, obj.peptide_list{i_list})
                % Extract information from the segment
                cur_ch = str2double(segment{5});        % charge
                cur_mz = str2double(segment{7})/cur_ch + CConstant.pmass;  % m/z
                MS2ScanI = str2double(segment{3});      % scan number
                
                % Find MS2 index
                mgf_name = erase(segment{2}, '.mgf');
                MS2_index = ms12DatasetIO.m_mapNameMS2Index(mgf_name);
                
                % Find the corresponding index of ms2
                tmp_idx = MS2_index(:,2) == MS2ScanI;
                if ~any(tmp_idx)
                    error('No matching MS2 scan found for scan number %d in %s', MS2ScanI, mgf_name); % No matching MS2 scan found
                end
                
                % Get MS1 scan number
                MS1Scan = MS2_index(tmp_idx, 1);
                
                % Get MS1 data
                MS1_index = ms12DatasetIO.m_mapNameMS1Index(mgf_name);
                MS1_peaks = ms12DatasetIO.m_mapNameMS1Peaks(mgf_name);
                
                ino = find(MS1_index(:,1) == MS1Scan);
                if isempty(ino)
                    error('No matching MS1 scan found for scan number %d in %s', MS1Scan, mgf_name); % No matching MS1 scan found
                end
                
                cur_rt = MS1_index(ino, 2); % retention time
                
                % Get peak list for this MS1
                first_peak_idx = [1; MS1_index(1:size(MS1_index), 3)];
                IX = first_peak_idx(ino):first_peak_idx(ino+1)-1;
                mz = MS1_peaks(IX, 1);
                inten = MS1_peaks(IX, 2);
                
                % Calculate tolerance
                if obj.ms1_tolerance.isppm
                    ptol = obj.ms1_tolerance.value * cur_mz * 1e-6;
                else
                    ptol = obj.ms1_tolerance.value;
                end
                
                % Find intensity within tolerance
                cur_inten = max(inten(mz >= cur_mz-ptol & mz <= cur_mz+ptol));
                if isempty(cur_inten)
                    break; % No intensity found within tolerance
                end
                
                % Calculate peptide mass
                lfMass = obj.getPeptideMass(obj.peptide_list{i_list});
                
                % Add to quantification
                pep_quant{i_list} = pep_quant{i_list}.appendOneSpecQuant(segment{2},cur_rt, ...
                    cur_inten,cur_mz,cur_ch,{segment{10}},lfMass,1); %#ok<CCAT1>
                break;
            end
        end
    end
    
    progress_printer.last_update();
end