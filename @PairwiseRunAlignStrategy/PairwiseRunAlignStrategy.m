classdef PairwiseRunAlignStrategy < CRunAlignStrategy
    % Pairwise alignment strategy with user-specified pairs

    properties
        m_pair_list
    end

    methods
        function obj = PairwiseRunAlignStrategy(pair_list)
            if nargin < 1
                pair_list = cell(0,2);
            end
            obj.m_pair_list = pair_list;
        end

        pairs = getAlignmentPairs(obj, raw_names)
    end
end
