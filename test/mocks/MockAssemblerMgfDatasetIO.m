classdef MockAssemblerMgfDatasetIO
    methods
        function [expPeaks, iCharge, precursorMZ] = read_oneSpec(~, ~, ~)
            expPeaks = [100, 200; 150, 100];
            iCharge = 2;
            precursorMZ = 500;
        end
    end
end
