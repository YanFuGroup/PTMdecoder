function drawGather(obj, pep_rtrange_map, dir_save, color_map, legend_map)
% Draw the XIC for gathered peptides using manually-checked rt range
% Input:
%   obj (CIMPGatherQuant)
%       Quantification aggregator instance
%   pep_rtrange_map (containers.Map)
%       map of [modified peptide _ charge _ raw file name] -> [rt_start, rt_end, check_label]
%   dir_save (1 x 1 char/string)
%       directory to save the figures
%   color_map (containers.Map or [])
%       color map (key: imp name, value: RGB 1x3)
%   legend_map (containers.Map or [])
%       legend map (key: imp name, value: display string)

% Check the input arguments
if nargin < 5
    legend_map = [];
end
if nargin < 4
    color_map = [];
end


% Do the same operation for gathered PSM for every raw file
keys_raw = obj.m_mapRawNames.keys;
for idx_keys = 1:obj.m_mapRawNames.Count
    idx_r = obj.m_mapRawNames(keys_raw{idx_keys});
    raw = obj.get_raw(idx_r);

    % Cluster the IMPs according to their masses
    group_idxs = cluster_imps_by_mass(raw.impMass,obj.m_ms1_tolerance);

    % Quantify the IMPs in each group
    for idx_g = 1:length(group_idxs)
        group_imp_name = raw.impNames(group_idxs{idx_g});
        group_ratio = raw.ratioMatrix(1:raw.length,group_idxs{idx_g});
        idxs_rt_inten = find(sum(group_ratio,2));
        group_ratio = group_ratio(idxs_rt_inten,:);
        group_rts = raw.curRts(idxs_rt_inten);
        group_inten = raw.curIntens(idxs_rt_inten);
        group_imp_mass = raw.impMass(group_idxs{idx_g});
        group_charge = raw.curCharge(idxs_rt_inten);
        [low_mz_bound, high_mz_bound, selected_charge, charge_group_idxs] = ...
            obj.get_mz_bound(group_imp_mass,group_charge);

        for idx_ch = 1:length(selected_charge)
            % Get retention time range for each IMP
            current_imp_rt_range = cell(length(group_imp_name),1);
            for idx_imp = 1:length(group_imp_name)
                generated_key = [group_imp_name{idx_imp},'_+', ...
                    num2str(selected_charge(idx_ch)), '_', keys_raw{idx_keys}];
                if pep_rtrange_map.isKey(generated_key)
                    current_imp_rt_range{idx_imp} = pep_rtrange_map(generated_key);
                end
            end
            if all(cellfun(@isempty,current_imp_rt_range))
                % All of the IMPs are removed in manual checking
                continue;
            end

            % Draw for this group
            obj.draw_each_group(keys_raw{idx_keys},...
                group_ratio(charge_group_idxs{idx_ch},:),...
                group_rts(charge_group_idxs{idx_ch},:),...
                group_inten(charge_group_idxs{idx_ch},:),...
                low_mz_bound(idx_ch),high_mz_bound(idx_ch), ...
                selected_charge(idx_ch),current_imp_rt_range,...
                group_imp_name, dir_save, color_map, legend_map);
        end
    end
end
end