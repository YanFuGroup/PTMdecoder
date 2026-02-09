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
                obj.blocks = CIMPQuantBlock.empty(0,1);
                return;
            end
            obj.blocks = blocks;
        end

        function obj = append_block(obj, block)
            % Append a CIMPQuantBlock to the report
            if isempty(block)
                return;
            end
            if isempty(obj.blocks)
                obj.blocks = block;
            else
                obj.blocks(end+1,1) = block;
            end
        end

        function obj = add_block(obj, protein_name_pos, records)
            % Convenience method to append a block using raw inputs
            block = CIMPQuantBlock(protein_name_pos, records);
            obj = obj.append_block(block);
        end
    end

    methods (Static)
        obj = empty(varargin)
    end

    methods
        
        pep_rtrange_map = build_pep_rtrange_map(obj)

        pep_prot_map = build_pep_prot_map(obj)
    end

end
