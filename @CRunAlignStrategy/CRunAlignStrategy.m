classdef (Abstract) CRunAlignStrategy < handle
    % Base class for run alignment strategies

    methods (Abstract)
        pairs = getAlignmentPairs(obj, raw_names)
    end
end
