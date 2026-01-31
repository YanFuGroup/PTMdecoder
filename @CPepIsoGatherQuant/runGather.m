function runGather(obj)
% Main entry point for summarizing the quantification of various IMPs of a peptide
% Output to a file
% Below each peptide, there are several lines starting with '*', representing the overall information of IMP.
% Each line starting with '*' (IMP line) is followed by several retention time lines ('@' starting lines).
%
% Input:
%   obj (CPepIsoGatherQuant)
%       Quantification aggregator instance

fout = fopen(obj.m_outputPath,'a');
if fout == -1
    error(['Cannot open the the report file ',obj.m_outputPath]);
end
has_written_protein = false;

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
        current_imp_name = obj.m_cstrIMPNames{idx_r}(group_idxs{idx_g});
        ratio_current = obj.m_ratioMatrix{idx_r}(:,group_idxs{idx_g});
        idxs_rt_inten = find(sum(ratio_current,2));
        ratio_current = ratio_current(idxs_rt_inten,:);
        current_rts = obj.m_curRts{idx_r}(idxs_rt_inten);
        current_inten = obj.m_curIntens{idx_r}(idxs_rt_inten);
        current_imp_mass = obj.m_IMPMass{idx_r}(group_idxs{idx_g});
        current_charge = obj.m_curCharge{idx_r}(idxs_rt_inten);
        [low_mz_bound, high_mz_bound, selected_charge, charge_group_idxs] = ...
            obj.get_mz_bound(current_imp_mass,current_charge);

        for idx_ch = 1:length(selected_charge)
            % group quant
            [has_nonzero_imp, imp_idx_nonzero, area_imp_final, rt_bound, idx_selected, ratio_each_XIC_peak] = ...
                obj.quant_each_group(keys_raw{idx_keys},...
                ratio_current(charge_group_idxs{idx_ch},:),...
                current_rts(charge_group_idxs{idx_ch},:),...
                current_inten(charge_group_idxs{idx_ch},:),...
                low_mz_bound(idx_ch),high_mz_bound(idx_ch),selected_charge(idx_ch));

            % Save to file
            % only write the non zero result
            if ~has_nonzero_imp
                % TODO?
                % Show that there is no quantification for this group
                continue;
            end
            if ~has_written_protein
                has_written_protein = true;
                write_protein_start_position_line(fout, obj.m_prot_names_pos);
            end
            for i_iso = 1:length(imp_idx_nonzero)
                % Write the line of peptide, charge, mass, area
                fprintf(fout, '*\t%s\t+%d\t%s\t%.4f\t%f\t%f\t%f\n', ...
                    current_imp_name{imp_idx_nonzero(i_iso)}, ...
                    selected_charge(idx_ch),...
                    keys_raw{idx_keys},...
                    mean([low_mz_bound(idx_ch),high_mz_bound(idx_ch)]),...
                    low_mz_bound(idx_ch),...
                    high_mz_bound(idx_ch),...
                    area_imp_final(i_iso,1));
                % Write the line of rt range, ratio, check label
                for i_peak = 1:size(rt_bound,2)
                    % If the ratio is 0, skip
                    if ratio_each_XIC_peak(i_iso,i_peak) == 0
                        continue;
                    end
                    if i_peak == idx_selected(i_iso)
                        check_label = 1;
                    else
                        check_label = 0;
                    end
                    fprintf(fout, '@\t%f\t%f\t%f\t%d\n',...
                        rt_bound(i_iso,i_peak).start,...
                        rt_bound(i_iso,i_peak).end, ...
                        ratio_each_XIC_peak(i_iso,i_peak), ...
                        check_label);
                end
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
