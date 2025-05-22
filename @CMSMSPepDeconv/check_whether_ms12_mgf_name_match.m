function check_whether_ms12_mgf_name_match(obj)
% Check whether the name of ms1/ms2 is matched with mgf

mgf_files = dir(fullfile(obj.m_specPath, '*.mgf'));
ms1_files = dir(fullfile(obj.m_specPath, '*.ms1'));

% Get the stem name of the files
mgf_names = {mgf_files.name};
ms1_names = {ms1_files.name};

% Remove the extension
mgf_names = cellfun(@(x) x(1:end-4), mgf_names, 'UniformOutput', false);
ms1_names = cellfun(@(x) x(1:end-4), ms1_names, 'UniformOutput', false);

% Check whether each .mgf file has a .ms1 file with the same name
for i = 1:length(mgf_names)
    if ~ismember(mgf_names{i}, ms1_names)
        error('File %s.mgf does not have a corresponding .ms1 file.', mgf_names{i});
    end
end
end