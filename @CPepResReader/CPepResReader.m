classdef CPepResReader
    % Read the result of peptide level
    % Usage:
    %   pep_reader = CPepResReader();
    %   pep_reader = pep_reader.read_from_pep_res_file(pep_res_path);
    
    properties
        m_pep_res;       % The overall peptide results, structured as:
        %   Map('mod_peptide_key':'charge',{}'dataset_name',{},'mean_mz',{},'lb_mz',{},'ub_mz',{},'rt_ranges',{},'ori_pep_line',{});
        m_prot_pep_res;  % The peptide results for each protein, structured as:
        %   cell('protein_name','mod_peptide_key');
    end

    methods
        function obj = CPepResReader()
            % Initialize the reader with an empty map
            obj.m_pep_res = containers.Map();
        end

        % Append one peptide result
        obj = append_one_pep(obj, key_mod_pep, charge_state, dataset_name, mean_mz, lb_mz, ub_mz, rt_ranges_temp, ori_pep_line);

        % Read from a peptide result file
        obj = read_from_pep_res_file(obj, pep_res_path);

        % Get retention time ranges for a modified peptide
        rt_ranges = get_rt_ranges(obj, mod_peptide);

        % Get the map of modified peptide to retention time ranges
        pep_rtrange_map = get_pep_rtrange_map(obj);
    end
end