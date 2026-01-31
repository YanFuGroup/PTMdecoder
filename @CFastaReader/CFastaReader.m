classdef CFastaReader
    % Read sequence library
    
    properties
        m_strFilePath;  % fasta file path
        m_regular_exp;  % regular expression to parse the protein name
    end
    
    methods
        function obj = CFastaReader(strFilePath, regular_exp)
            % Input:
            %   strFilePath (1 x 1 char/string)
            %       FASTA file path
            %   regular_exp (1 x 1 char/string)
            %       regular expression to parse the protein name
            obj.m_strFilePath = strFilePath;
            obj.m_regular_exp = regular_exp;
        end
        
        % Read sequence from file to m_mapProt
        mapProt = read(obj);
    end
end
