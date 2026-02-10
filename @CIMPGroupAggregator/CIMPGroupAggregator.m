classdef CIMPGroupAggregator < handle
    % Aggregate IMP groups across raw/charge and invoke a callback

    properties (Access=private)
        m_ms1_tolerance
    end

    methods
        % Constructor for CIMPGroupAggregator
        function obj = CIMPGroupAggregator(ms1_tolerance)
            obj.m_ms1_tolerance = ms1_tolerance;
        end

        % Aggregate IMP groups across raw files and invoke callback
        state = aggregate(obj, raw_names, raw_ident_stores, pep_rtrange_map, handle_group_charge, state);
    end

    methods (Access=private)
        % Calculate m/z bounds for IMP groups
        [low_mz_bound, high_mz_bound, selected_charge, charge_group_idxs] = ...
            getMzBound(obj, current_imp_mass, current_charge);
        
        % Build retention time ranges for IMPs
        current_imp_rt_range = buildImpRtRanges(obj, group_imp_name, selected_charge, raw_name, pep_rtrange_map);
        
        % Cluster the IMPs according to their masses.
        idxs_res = cluster_by_mass(obj, imp_masses);
    end
end
