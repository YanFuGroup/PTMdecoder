classdef CIMPQuantStats
    % Statistics helpers for IMP quantification routines

    methods (Static)
        % Record or flush rt_sorted lengths for offline analysis
        rt_sorted_stats(action, value)

        % Record or summarize quantification group outcomes
        stats = quant_group_stats(action, value)

        % Write an INFO summary for initial quantification group outcomes
        log_quant_group_summary(stage_name)
    end
end
