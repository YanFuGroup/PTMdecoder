classdef CXICAlignAnchorSelector
    % Select unmodified peptide anchors from FDR filtered results

    methods
        anchors = selectAnchors(obj, fdr_filtered_result_path, ms12DatasetIO, options)
    end
end
