function processedPeaks = processPeaks(expPeaks, alpha)
% processPeaks - Normalize and remove low relative intensity peaks

processedPeaks = [expPeaks,expPeaks(:,2)];
if isempty(processedPeaks)
    return;
end

dFactor = max(processedPeaks(:,2));
processedPeaks(:,2) = processedPeaks(:,2) / (dFactor + eps);
processedPeaks(processedPeaks(:,2) < alpha,:) = [];
end
