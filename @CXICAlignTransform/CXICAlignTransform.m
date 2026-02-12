classdef CXICAlignTransform
    % Fit and apply RT alignment transforms

    methods (Static)
        model = fit(ref_rts, target_rts, options)

        rt_out = apply(model, rt_in)
    end
end
