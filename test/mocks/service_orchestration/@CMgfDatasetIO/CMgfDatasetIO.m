classdef CMgfDatasetIO
    methods
        function obj = CMgfDatasetIO(varargin)
        end


        function [expPeaks, iCharge, precursorMZ] = read_oneSpec(obj, datasetName, specName) %#ok<INUSD>
            expPeaks = [100, 200; 150, 100];
            iCharge = 2;
            precursorMZ = 500;
        end
    end
end
