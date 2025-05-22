function obj = append_one_spectrum(obj, dataset_name, spectrum_name)
% Add one spectrum at the end of the spectrum list
% Input:
%   spectrum_name
%       the spectrum name

obj.m_peps_specs_forms(end).spectrum_list(end+1).dataset_name = dataset_name;
obj.m_peps_specs_forms(end).spectrum_list(end).spectrum_name = spectrum_name;
obj.m_peps_specs_forms(end).spectrum_list(end).peptidoform_list_str = cell(50,1);
obj.m_peps_specs_forms(end).spectrum_list(end).peptidoform_list_abun = zeros(50,1);
obj.m_peps_specs_forms(end).spectrum_list(end).peptidoform_num = 0;

end

