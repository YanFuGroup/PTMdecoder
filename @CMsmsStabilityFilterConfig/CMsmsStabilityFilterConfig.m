classdef CMsmsStabilityFilterConfig
    % Config parser for peptide-level MSMS stability filter.

    methods (Static)
        filter_cfg = fromParamMap(task_param_map, err_id_prefix)
    end
end
