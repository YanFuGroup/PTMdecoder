function obj = append_one_peptide(obj, peptide_sequence)
% Add one peptide at the end of the peptide list
% Input:
%   peptide_sequence
%       the peptide sequence

obj.m_peps_specs_forms(end+1).peptide_sequence = peptide_sequence;
obj.m_peps_specs_forms(end).spectrum_list = struct('dataset_name',{},'spectrum_name',{}, ...
    'peptidoform_list_str',{},'peptidoform_list_abun',{},'peptidoform_num',{});

end

