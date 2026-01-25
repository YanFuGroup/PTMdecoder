classdef CPepProtService
    % Service for Protein and Peptide logic
    
    properties
        m_mapProt;      % Map of protein name to sequence
        m_pep_prot_mapper; % CPeptideProteinMap object for peptide to protein mapping
    end
    
    methods
        function obj = CPepProtService(fastaFilePath, regularExp, filtered_res_file_path)
            % Internalizes the reader creation logic (Composition)
            % This makes the service easier to use, as the caller only needs to provide file paths.
            reader = CFastaReader(fastaFilePath, regularExp);
            obj.m_mapProt = reader.read();
            
            obj.m_pep_prot_mapper = CPeptideProteinMap(filtered_res_file_path);
        end

        % Check whether it is the N-terminus or C-terminus of the protein
        [isProtN,isProtC] = getWhetherProtNC(obj,strSeq);

        % find the protein of a peptide
        cell_prot_name_pos = get_protein_name_pos(obj,pepSeq);
    end
end
