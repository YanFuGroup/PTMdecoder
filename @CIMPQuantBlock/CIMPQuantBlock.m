classdef CIMPQuantBlock
    % Data model for a protein block in an IMP quantification report

    properties
        protein_name_pos
        records
    end

    methods
        function obj = CIMPQuantBlock(protein_name_pos, records)
            % Constructor for CIMPQuantBlock
            % input:
            %   protein_name_pos (N x 2 cell)
            %       protein names and their start positions
            %   records (CIMPQuantRecord array)
            %       IMP quantification records for this block
            if nargin == 0
                return;
            end
            obj.protein_name_pos = protein_name_pos;
            obj.records = records;
        end
    end

    methods (Static)
        function obj = empty(varargin)
            obj = CIMPQuantBlock;
            obj = obj([]);
        end
    end
end
