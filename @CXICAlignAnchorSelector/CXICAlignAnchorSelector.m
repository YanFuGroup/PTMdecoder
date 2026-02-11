classdef CXICAlignAnchorSelector
    % Select unmodified peptide anchors from FDR filtered results

    methods
        anchors = selectAnchors(obj, fdr_file_path, ms12DatasetIO, options)
    end
end
