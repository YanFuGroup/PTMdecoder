function [low_mz_bound, high_mz_bound, selected_charge, charge_group_idxs] = ...
    get_mz_bound(obj,current_iso_mass,current_charge)
% Get the m/z bound of ms1 peak
% input:    current_iso_mass
%               current neutral peptide IMP mass
%           current_charge
%               current precursor charge
% output:   low_bound
%               low bound of m/z ms1 peak
%           high_bound
%               high bound of m/z ms1 peak
%           selected_charge
%               current precursor charge
%           (every output is in column matrix form, each row corresponds to a charge state)
%           group_idxs
%               indexes of each charge, in cell form

% Record the high bound and the low bound of observed precursor m/z.
% Return different mz bounds for different charge state
[selected_charge, ~, group_belong] = unique(current_charge);
high_mz_bound = zeros(length(selected_charge),1);
low_mz_bound = zeros(length(selected_charge),1);
charge_group_idxs = cell(length(selected_charge),1);
for idx_ch = 1:length(selected_charge)
    charge_group_idxs{idx_ch} = find(group_belong==idx_ch);
    temp_charge = selected_charge(idx_ch);
    high_mz_bound(idx_ch) = max((current_iso_mass+temp_charge*CConstant.pmass)./temp_charge);
    low_mz_bound(idx_ch) = min((current_iso_mass+temp_charge*CConstant.pmass)./temp_charge);
end

% % Record the high bound and the low bound of observed precursor m/z.
% % Return different mz bounds for different charge state
% [selected_charge, ~, group_belong] = unique(current_charge);
% high_mz_bound = zeros(length(selected_charge),1);
% low_mz_bound = zeros(length(selected_charge),1);
% charge_group_idxs = cell(length(selected_charge),1);
% for idx_ch = 1:length(selected_charge)
%     charge_group_idxs{idx_ch} = find(group_belong==idx_ch);
%     temp_charge = selected_charge(idx_ch);
%     high_mz_bound(idx_ch) = max(current_mz(temp_charge == current_charge));
%     low_mz_bound(idx_ch) = min(current_mz(temp_charge == current_charge));
% end

% Add tolerance to high bound, and minus tolerance to low bound.
if obj.m_ms1_tolerance.isppm
    high_mz_bound = high_mz_bound .* (1+obj.m_ms1_tolerance.value*1e-6);
    low_mz_bound = low_mz_bound .* (1-obj.m_ms1_tolerance.value*1e-6);
else
    high_mz_bound = high_mz_bound + obj.m_ms1_tolerance.value;
    low_mz_bound = low_mz_bound - obj.m_ms1_tolerance.value;
end
end