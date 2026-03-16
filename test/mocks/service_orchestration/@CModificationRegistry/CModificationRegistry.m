classdef CModificationRegistry
    methods (Static)
        function [fixedModNameMass, variableModNameMass] = fromConfig(~)
            fixedModNameMass = containers.Map('KeyType', 'char', 'ValueType', 'double');
            variableModNameMass = containers.Map('KeyType', 'char', 'ValueType', 'double');
        end
    end
end
