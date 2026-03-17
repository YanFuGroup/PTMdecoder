classdef CMS2QuantSolver
    methods (Static)
        function noise_model = estimateDatasetNoiseModel(noise_model_fit_inputs) %#ok<INUSD>
            noise_model = struct('sigma_base', 1.0, 'gamma', 0.1, 'tau_floor', 5.0);
        end


        function stability_diag = estimateStability(~, ~, ~, ~, base_abundance, ~, ~, stability_options)
            reported_imp_indices = find(base_abundance > 0);
            support_frequency = 0.8 * ones(numel(reported_imp_indices), 1);
            abundance_mad = 0.05 * ones(numel(reported_imp_indices), 1);
            stability_diag = struct( ...
                'jaccard_stability', 0.9, ...
                'support_frequency', support_frequency, ...
                'abundance_mad', abundance_mad, ...
                'reported_imp_indices', reported_imp_indices, ...
                'num_successful_resamples', stability_options.n_resamples);
        end
    end
end
