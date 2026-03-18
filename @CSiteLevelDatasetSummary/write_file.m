function write_file(obj)
% Write site-by-dataset matrix to file.
% Inputs:
%   obj (CSiteLevelDatasetSummary)
%       dataset-level site summarizer instance

[output_dir, ~, ~] = fileparts(obj.m_output_site_dataset_matrix_path);
CPathResolver.ensureDir(output_dir);

fout = fopen(obj.m_output_site_dataset_matrix_path, 'w');
if fout < 0
    error(['The file "', obj.m_output_site_dataset_matrix_path, '" is not available.']);
end

% Header: Site + dataset columns
fprintf(fout, 'Site');
for idx_dataset = 1:numel(obj.m_dataset_names)
    fprintf(fout, '\t%s', obj.m_dataset_names{idx_dataset});
end
fprintf(fout, '\n');

% Matrix rows
for idx_site = 1:numel(obj.m_site_names)
    site_name = obj.m_site_names{idx_site};
    fprintf(fout, '%s', site_name);

    dataset_sum_map = obj.m_site_dataset_sum(site_name);
    for idx_dataset = 1:numel(obj.m_dataset_names)
        dataset_name = obj.m_dataset_names{idx_dataset};
        site_dataset_value = 0;
        if isKey(dataset_sum_map, dataset_name)
            site_dataset_value = dataset_sum_map(dataset_name);
        end
        fprintf(fout, '\t%.15g', site_dataset_value);
    end
    fprintf(fout, '\n');
end

fclose(fout);

CLogger.info(['[CSiteLevelDatasetSummary:write_file] Matrix written. ', ...
    'output=%s, site_count=%d, dataset_count=%d.'], ...
    obj.m_output_site_dataset_matrix_path, numel(obj.m_site_names), numel(obj.m_dataset_names));
end
