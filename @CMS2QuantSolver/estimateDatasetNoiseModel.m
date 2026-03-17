function noise_model = estimateDatasetNoiseModel(noise_model_fit_inputs)
% estimateDatasetNoiseModel - Estimate global noise parameters for MS2 matched peaks
% Learn dataset-level global noise parameters using pooled MLE and MAD robust 
% estimation to create a shared noise model for subsequent stability perturbations.
% Inputs:
%   noise_model_fit_inputs (1 x N struct)
%       Struct array from first phase baseline solves. Required fields for each element: 
%       - filteredOutExpPeakCount (1 x 1 double):
%           Number of background peaks in normalized-intensity space.
%       - filteredOutExpPeakSqSum (1 x 1 double):
%           Sum of squared normalized intensities of background peaks.
%       - matchedExpPeaks (M x 3 double):
%           Matrix of matched peaks (column 2 is normalized intensity).
%       - fittedMatchedPeakIntensities (M x 1 double):
%           Reconstructed intensities (yhat).
% Outputs:
%   noise_model (1 x 1 struct)
%       Struct with global parameters: 
%       - sigma_base (1 x 1 double):
%           Baseline noise standard deviation. Set to eps if no background.
%       - gamma (1 x 1 double):
%           Multiplicative scaling factor. Set to 0 if no high signals.
%       - tau_floor (1 x 1 double):
%           Minimum floor threshold (5 * sigma_base).

if nargin < 1 || isempty(noise_model_fit_inputs)
    CLogger.error('CMS2QuantSolver:estimateDatasetNoiseModel:EmptyInput', 'The input noise_model_fit_inputs is empty or undefined.');
end

% 1. Estimate sigma_base (Pooled MLE)
totalBackgroundCount = sum([noise_model_fit_inputs.filteredOutExpPeakCount]);
totalBackgroundSqSum = sum([noise_model_fit_inputs.filteredOutExpPeakSqSum]);

if totalBackgroundCount == 0
    CLogger.warn('CMS2QuantSolver:estimateDatasetNoiseModel:NoBackground', 'Total background count is 0, setting sigma_base to eps.');
    sigma_base = eps;
else
    sigma_base = sqrt(totalBackgroundSqSum / totalBackgroundCount);
end

% 2. Threshold setting
tau_floor = max(5 * sigma_base, eps);

% 3. Collect residuals for gamma estimation
all_q_cell = cell(1, numel(noise_model_fit_inputs));

for i = 1:numel(noise_model_fit_inputs)
    inputs_i = noise_model_fit_inputs(i);

    yhat = inputs_i.fittedMatchedPeakIntensities;
    if isempty(yhat)
        continue;
    end
    
    if size(inputs_i.matchedExpPeaks, 2) < 2
        CLogger.error('CMS2QuantSolver:estimateDatasetNoiseModel:InvalidMatchedPeaks', ...
            'matchedExpPeaks must contain normalized intensity in column 2.');
    end
    y = inputs_i.matchedExpPeaks(:, 2);
    
    % Filter high signal peaks
    maskHigh = yhat > tau_floor;
    if ~any(maskHigh)
        continue;
    end
    
    yHigh = y(maskHigh);
    yhatHigh = yhat(maskHigh);
    
    % Calculate residuals
    r_i = yHigh - yhatHigh;
    
    % Calculate relative residuals q_i
    % The maskHigh filter already ensures yhatHigh > tau_floor, so max() is redundant.
    q_i = r_i ./ yhatHigh;

    all_q_cell{i} = q_i(:);
end

all_q = vertcat(all_q_cell{:});

% 4. Estimate gamma (Robust MAD)
if isempty(all_q)
    CLogger.warn('CMS2QuantSolver:estimateDatasetNoiseModel:NoHighSignals', 'No high signal peaks found across dataset. Setting gamma to 0.');
    gamma = 0;
else
    % Robust scale estimate for high-signal relative residuals.
    % Under a zero-centered residual assumption, use the sample median of q.
    mad_q = median(abs(all_q));
    gamma = 1.4826 * mad_q;
end

% 5. Build output noise model
noise_model.sigma_base = sigma_base;
noise_model.gamma = gamma;
noise_model.tau_floor = tau_floor;

end
