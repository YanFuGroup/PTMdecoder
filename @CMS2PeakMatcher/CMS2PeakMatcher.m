classdef CMS2PeakMatcher
    % Static matching and post-match filtering for MS2 peaks

    methods (Static)
        matchedExpPeaks = match(expPeaks, vNonRedunTheoryIonMz, ms2_tolerance)

        processedPeaks = processPeaks(expPeaks, alpha)

        [massArrangement, vNonRedunTheoryIonMz] = postMatchFilterPeptidoforms(matchedExpPeaks, massArrangement, vNonRedunTheoryIonMz)
    end
end
