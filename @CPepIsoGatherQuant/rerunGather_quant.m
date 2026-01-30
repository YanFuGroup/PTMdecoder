function rerunGather_quant(obj,pep_rtrange_map)
% Re-quantification for gathered peptides using manually-checked rt range
% Input:
%   obj (CPepIsoGatherQuant)
%       Quantification aggregator instance
%   pep_rtrange_map (containers.Map)
%       map of [modified peptide _ charge _ raw file name] -> [rt_start, rt_end, check_label]

fout = fopen(obj.m_outputPath,'a');
if fout == -1
    error(['Cannot open the the report file ',obj.m_outputPath]);
end
is_protein_writen = false;

% Delete useless raws which are empty or blank.
obj = obj.delUselessRaws();

% Do the same operation for gathered PSM for every raw file
keys_raw = obj.m_mapRawNames.keys;
for idx_keys = 1:obj.m_mapRawNames.Count
    idx_r = obj.m_mapRawNames(keys_raw{idx_keys});

    % Cluster the IMPs according to their masses
    group_idxs = clustering_IMPs(obj.m_IMPMass{idx_r},obj.m_ms1_tolerance);

    % Quantify the IMPs in each group
    for idx_g = 1:length(group_idxs)
        current_iso_name = obj.m_cstrIMPNames{idx_r}(group_idxs{idx_g});
        current_ratioMatrix = obj.m_ratioMatrix{idx_r}(:,group_idxs{idx_g});
        idxs_rt_inten = find(sum(current_ratioMatrix,2));
        current_ratioMatrix = current_ratioMatrix(idxs_rt_inten,:);
        current_rts = obj.m_curRts{idx_r}(idxs_rt_inten);
        current_inten = obj.m_curIntens{idx_r}(idxs_rt_inten);
        current_iso_mass = obj.m_IMPMass{idx_r}(group_idxs{idx_g});
        current_charge = obj.m_curCharge{idx_r}(idxs_rt_inten);
        [low_mz_bound, high_mz_bound, selected_charge, charge_group_idxs] = ...
            obj.get_mz_bound(current_iso_mass,current_charge);

        for idx_ch = 1:length(selected_charge)
            % Get retention time range for each IMP
            current_iso_rt_range = cell(length(current_iso_name),1);
            for idx_imp = 1:length(current_iso_name)
                generated_key = [current_iso_name{idx_imp},'_+', ...
                    num2str(selected_charge(idx_ch)), '_', keys_raw{idx_keys}];
                if pep_rtrange_map.isKey(generated_key)
                    current_iso_rt_range{idx_imp} = pep_rtrange_map(generated_key);
                end
            end
            if all(cellfun(@isempty,current_iso_rt_range))
                % All of the IMPs are removed in manual checking
                continue;
            end

            % group quant
            [bhave_non_zeros, idx_cur, area, rt_bound, max_label, ratio_each_XIC_peak] = ...
                obj.requant_each_group(keys_raw{idx_keys},...
                current_ratioMatrix(charge_group_idxs{idx_ch},:),...
                current_rts(charge_group_idxs{idx_ch},:),...
                current_inten(charge_group_idxs{idx_ch},:),...
                low_mz_bound(idx_ch),high_mz_bound(idx_ch), ...
                selected_charge(idx_ch),current_iso_rt_range);

            % Save to file
            % only write the non zero result
            if ~bhave_non_zeros
                continue;
            end
            if ~is_protein_writen
                is_protein_writen = true;
                write_protein_start_position_line(fout, obj.m_prot_names_pos);
            end
            for i_iso = 1:length(idx_cur)
                fprintf(fout, '*\t%s\t+%d\t%s\t%.4f\t%f\t%f\t%f\n', ...
                    current_iso_name{idx_cur(i_iso)}, ...
                    selected_charge(idx_ch),...
                    keys_raw{idx_keys},...
                    mean([low_mz_bound(idx_ch),high_mz_bound(idx_ch)]),...
                    low_mz_bound(idx_ch),...
                    high_mz_bound(idx_ch),...
                    area(i_iso,1));
                fprintf(fout, '@\t%f\t%f\t%f\t%d\n', ...
                    rt_bound(i_iso).start, ...
                    rt_bound(i_iso).end, ...
                    ratio_each_XIC_peak(i_iso), ...
                    max_label(i_iso));
            end
        end
    end
end
fclose(fout);
end

% Cluster the IMPs according to their masses
function idxs_res = clustering_IMPs(IMP_masses,ms1_tolerance)
% Input:
%   IMP_masses (1 x K double) Da
%       Masses of IMPs
%   ms1_tolerance (struct)
%       Tolerance of MS1 (fields: isppm, value)
% Output:
%   idxs_res (1 x G cell)
%       Indices of each group, in cell form
[m_val,m_inx] = sort(IMP_masses);

% Initialize variables
idxs_res = {};
currentCluster = [];

% Iterate through the numbers
for i = 1:length(m_val)-1
    currentNumber = m_val(i);
    nextNumber = m_val(i+1);
    % Add the current number to the current cluster
    currentCluster = [currentCluster, m_inx(i)]; %#ok<AGROW>

    % Mass tolerance
    if ms1_tolerance.isppm
        tol = (ms1_tolerance.value * currentNumber)/1e6;
    else
        tol = ms1_tolerance.value;
    end
    
    % Check if the next mass is beyond the threshold
    if abs(nextNumber - currentNumber) > tol
        % Add the current cluster to the list of clusters
        idxs_res = [idxs_res, {currentCluster}]; %#ok<AGROW> 
        
        % Reset the current cluster
        currentCluster = [];
    end
end

% Add the last number to the current cluster
currentCluster = [currentCluster, m_inx(end)];

% Add the last cluster to the list of clusters
idxs_res = [idxs_res, {currentCluster}];  

end



% Write protein start position line
function write_protein_start_position_line(fid, prot_names_pos)
% Input:
%   fid (1 x 1 double/int)
%       File identifier
%   prot_names_pos (P x 2 cell)
%       Protein name and start position pairs
for idx_np = 1:size(prot_names_pos,1)
    fprintf(fid, '%s,%d;', prot_names_pos{idx_np,1},...
        prot_names_pos{idx_np,2});
end
fprintf(fid,'\n');
end