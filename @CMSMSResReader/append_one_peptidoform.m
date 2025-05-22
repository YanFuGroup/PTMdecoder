function obj = append_one_peptidoform(obj, peptidoform_str, relative_abundance)
% Add one peptidoform at the end of the peptidoform list
% Input:
%   peptidoform_str
%       the peptidoform string
%   relative_abundance
%       the relative abundance of the peptidoform

obj.m_peps_specs_forms(end).spectrum_list(end).peptidoform_num = ...
    obj.m_peps_specs_forms(end).spectrum_list(end).peptidoform_num + 1;
len = obj.m_peps_specs_forms(end).spectrum_list(end).peptidoform_num;
% If the length is larger than the buffer, extend the buffer
if len > length(obj.m_peps_specs_forms(end).spectrum_list(end).peptidoform_list_abun)
    obj.m_peps_specs_forms(end).spectrum_list(end).peptidoform_list_str{end+50} = '';
    obj.m_peps_specs_forms(end).spectrum_list(end).peptidoform_list_abun(end+50) = 0;
end

obj.m_peps_specs_forms(end).spectrum_list(end).peptidoform_list_str{len} = peptidoform_str;
obj.m_peps_specs_forms(end).spectrum_list(end).peptidoform_list_abun(len) = relative_abundance;

end

