classdef CMSMSResReader
    % Read the result of msms level
    % Usage:
    %   msms_reader = CMSMSResReader();
    %   msms_reader = msms_reader.read_from_msms_res_file(msms_res_path);
    %   disp(msms_reader.m_peps_specs_forms(i_pep).peptide_sequence);
    %   disp(msms_reader.m_peps_specs_forms(i_pep).spectrum_list(i_spec).dataset_name);
    %   disp(msms_reader.m_peps_specs_forms(i_pep).spectrum_list(i_spec).peptidoform_list_str{i_pls});
    %   disp(msms_reader.m_peps_specs_forms(i_pep).spectrum_list(i_spec).peptidoform_list_abun(i_pla));
    
    properties(Access = public)
        % The hierarchical structure of the overall MSMS results: 
        %   struct('peptide_sequence',{},'spectrum_list',{})
        %   - peptide_sequence: the sequence of the peptide
        %   - spectrum_list: the struct of the spectra of the peptide
        %       struct('dataset_name',{},'spectrum_name',{},'peptidoform_list_str',{},'peptidoform_list_abun',{},'peptidoform_num',{})
        %       - dataset_name: the name of the dataset
        %       - spectrum_name: the name of the spectrum
        %       - peptidoform_list_str: the list of the peptidoforms in the spectrum, cell array (peptidoform_num x 1)
        %       - peptidoform_list_abun: the list of the relative abundance of the peptidoforms in the spectrum, double array (peptidoform_num x 1)
        %       - peptidoform_num: the number of the peptidoforms in the spectrum
        m_peps_specs_forms;
    end
    
    methods
        function obj = CMSMSResReader()
            obj.m_peps_specs_forms = struct('peptide_sequence',{},'spectrum_list',{});
        end

        % Add one peptide
        obj = append_one_peptide(obj, peptide_sequence);

        % Add one spectrum
        obj = append_one_spectrum(obj, dataset_name, spectrum_name);

        % Add one peptidoform
        obj = append_one_peptidoform(obj, peptidoform_str, relative_abundance);

        % Read from a msms result file
        obj = read_from_msms_res_file(obj, msms_res_path);

    end
end