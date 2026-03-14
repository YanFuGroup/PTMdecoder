function fitted = computeFittedMatchedPeakIntensities(vNonRedunTheoryIonMz, matchedExpPeaks, abundance, frageff, ionTypePosCharge)
% Reconstructs the fitted matched-peak intensities (yhat_i).
% 
% Purpose:
%   Distinguishes between a primary path (using site-specific fragment efficiency
%   if both `frageff` and `ionTypePosCharge` are available) and a fallback path
%   (global scaling regression) when they are missing.
% 
% Note on ion metadata query:
%   This implementation continues using the legacy 3-column query
%   `[type, position, charge]` for matching `frageff`. A full 4-column
%   `ionGroupInfo` with `classIndex` is deferred to a future refactoring
%   for backwards compatibility.
%
% Input:
%   vNonRedunTheoryIonMz (L x T double)
%       Non-redundant theoretical ion table. Columns 7:end are indicator 
%       values for peptidoforms. (Cols 2:4 are [type, position, charge]).
%   matchedExpPeaks (K x 3 double)
%       Matrix representing experimental peaks matched to theoretical ions.
%       Columns: [ion_index, normalized_intensity, raw_intensity].
%   abundance (M x 1 double)
%       Estimated relative abundance vector for each peptidoform.
%   frageff (Q x 1 double, optional)
%       Estimated fragmentation efficiency coefficients.
%   ionTypePosCharge (Q x 3 double, optional)
%       Unique ion grouping keys formatted as `[type, position, charge]`.
%
% Output:
%   fitted (K x 1 double)
%       Reconstructed intensities corresponding to each matchedExpPeak row.

if isempty(matchedExpPeaks)
    fitted = [];
    return;
end

if nargin < 4
    frageff = [];
end
if nargin < 5
    ionTypePosCharge = [];
end

% 1. Extract matched peak information
peakIonIndices = matchedExpPeaks(:, 1);
expIntensities = matchedExpPeaks(:, 2);

% 2. Compute base structural intensities
% Columns 7:end match the M possible peptidoforms
impIndicators = vNonRedunTheoryIonMz(peakIonIndices, 7:end);
baseIntensities = impIndicators * abundance(:);

% Check if primary path information is fully available
hasPrimaryInfo = ~isempty(frageff) && ~isempty(ionTypePosCharge) && ...
                    (length(frageff) == size(ionTypePosCharge, 1));
                    
fitted = zeros(size(expIntensities));

if hasPrimaryInfo
    % === Primary Path ===
    % Extract matched ions metadata: [type, position, charge] (Cols 2:4)
    matchedIonMetadata = vNonRedunTheoryIonMz(peakIonIndices, 2:4);
    
    % Map each matched peak to its fragment efficiency coefficient
    [validIdx, frageffIdx] = ismember(matchedIonMetadata, ionTypePosCharge, 'rows');
    
    fittedShape = zeros(size(peakIonIndices));
    
    % Shape = Structural Base * Site-specific Efficiency
    fittedShape(validIdx) = baseIntensities(validIdx) .* frageff(frageffIdx(validIdx));
    
    % Align the locally-weighted shape back to experimental scale 
    % via non-negative least squares regression: c = max(0, shape \ exp).
    if any(fittedShape > 0)
        scalingFactor = max(0, fittedShape \ expIntensities);
        fitted = fittedShape * scalingFactor;
    end
else
    % === Fallback Path ===
    % Directly uses global scaling factor
    if any(baseIntensities > 0)
        scalingFactor = max(0, baseIntensities \ expIntensities);
        fitted = baseIntensities * scalingFactor;
    end
end
end
