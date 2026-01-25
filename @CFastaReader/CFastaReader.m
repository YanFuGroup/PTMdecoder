classdef CFastaReader
    % Read sequence library
    
    properties
        m_strFilePath;  % fasta file path
        m_regular_exp;  % regular expression to parse the protein name
    end
    
    methods
        function obj = CFastaReader(strFilePath, regular_exp)
            obj.m_strFilePath = strFilePath;
            obj.m_regular_exp = regular_exp;
        end
        
        % Read sequence from file to m_mapProt
        mapProt = read(obj);
    end
end
