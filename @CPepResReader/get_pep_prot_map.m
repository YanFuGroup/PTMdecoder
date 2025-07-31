function pep_prot_map = get_pep_prot_map(obj)
% Get the peptide-protein mapping
% Output:
%   pep_prot_map - A map where keys are modified peptide strings and values are cell arrays of protein names and sites

pep_prot_map = containers.Map();
for i = 1:size(obj.m_prot_pep_res, 1)
    protein_name = obj.m_prot_pep_res{i, 1};
    peptides = obj.m_prot_pep_res{i, 2};
    
    for j = 1:length(peptides)
        mod_peptide = peptides{j};
        pep_prot_map(mod_peptide) = protein_name;
    end
end