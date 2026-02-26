function rawIdentManager = buildRawIdentManagerFromSpectrumList(obj, spectrum_list)
% Build a raw identification manager from a spectrum list
% Input:
%   obj (CMSMSPepDeconv)
%       Processor instance
%   spectrum_list (struct array)
%       spectrum entries with fields:
%       - dataset_name
%       - spectrum_name
%       - peptidoform_list_str
%       - peptidoform_list_abun
% Output:
%   rawIdentManager (CIMPRawIdentManager)
%       per-raw identification store manager

deps = struct( ...
    'getProfilesFunc', @(dataset_name, spectrum_name) obj.getProfiles(dataset_name, spectrum_name), ...
    'fixedModNameMass', {obj.m_fixedModNameMass}, ...
    'variableModNameMass', {obj.m_variableModNameMass});

rawIdentManager = CPeptideRawIdentAssembler.buildFromSpectrumList(spectrum_list, deps);
end
