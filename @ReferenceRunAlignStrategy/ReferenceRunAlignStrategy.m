classdef ReferenceRunAlignStrategy < CRunAlignStrategy
    % Reference-run alignment strategy

    properties
        m_reference_raw
    end

    methods
        function obj = ReferenceRunAlignStrategy(reference_raw)
            if nargin < 1
                reference_raw = '';
            end
            obj.m_reference_raw = reference_raw;
        end

        pairs = getAlignmentPairs(obj, raw_names)
    end
end
