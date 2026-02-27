classdef CPeptideAlignRequantServiceConfig
    % Config builder for CPeptideAlignRequantService.

    methods (Static)
        cfg = fromParamMap(task_param_map)
    end

    methods (Static, Access = private)
        cfg = finalize(cfg)
    end
end
