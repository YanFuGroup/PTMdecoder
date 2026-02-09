classdef CIMPQuantReport
    % Data model for an IMP quantification report file

    properties
        blocks
    end

    methods
        function obj = CIMPQuantReport(blocks)
            % Constructor for CIMPQuantReport
            % input:
            %   blocks (CIMPQuantBlock array)
            %       protein blocks in the report
            if nargin == 0
                return;
            end
            obj.blocks = blocks;
        end
    end

    methods (Static)
        obj = empty(varargin)
    end

    methods
        write(obj, path)

        pep_rtrange_map = build_pep_rtrange_map(obj)

        pep_prot_map = build_pep_prot_map(obj)
    end

end
