classdef CFastaIO
    % Read sequence library
    
    properties
        m_strFilePath;  % fasta file path
        m_regular_exp;  % regular expression to parse the protein name
        m_mapProt;      % Map of protein name to sequence
        m_pep_prot_mapper; % CPeptideProteinMap object for peptide to protein mapping
        % m_strOneLine;   % Represent the entire fasta file as a single string
        % m_strProtName;  % the name of all proteins, using regular expression '>([^ ,]*)'
        % m_strProtPos;   % the position of all proteins in single line (m_strProtSeq)
    end
    
    methods
        function obj = CFastaIO(strFilePath, regular_exp, filtered_res_file_path)
            obj.m_strFilePath = strFilePath;
            obj.m_regular_exp = regular_exp;
            obj.m_pep_prot_mapper = CPeptideProteinMap(filtered_res_file_path);
            obj = Init(obj);
        end
        
        % Initialize, read sequence from file to m_strProtSeq
        obj = Init(obj);

        % Check whether it is the N-terminus or C-terminus of the protein
        [isProtN,isProtC] = getWhetherProtNC(obj,strSeq);

        % find the protein of a peptide
        cell_prot_name_pos = get_protein_name_pos(obj,pepSeq);
    end
end

