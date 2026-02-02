function proteinCell = get_proteins(obj, peptide)
    % get_proteins takes a peptide string as input and returns a cell array of corresponding protein names
    % Input:
    %   obj (CPeptideProteinMap)
    %       Peptide->protein mapping instance
    %   peptide (1 x 1 char/string)
    %       Peptide sequence
    % Output:
    %   proteinCell (1 x N cell)
    %       Cell array of protein names, {} if not found
    
    if ischar(peptide)
        pepKey = peptide;
    elseif isstring(peptide)
        pepKey = char(peptide);
    else
        proteinCell = {};
        return;
    end
    
    if isKey(obj.m_map, pepKey)
        proteinCell = obj.m_map(pepKey);
    else
        proteinCell = {};
    end
end