function compress(obj)
% Trim buffers and remove empty entries
% Organize the data structure: delete the empty peptide, spectrum, and peptidoform
% Also trim the buffers

% Iterate backwards to safely delete
for i = length(obj.Peptides):-1:1
    if isempty(obj.Peptides(i).spectrum_list)
        obj.Peptides(i) = [];
    else
        for j = length(obj.Peptides(i).spectrum_list):-1:1
            num = obj.Peptides(i).spectrum_list(j).peptidoform_num;
            % Check for 0 or empty (uninitialized buffer slots)
            if isempty(num) || num == 0
                obj.Peptides(i).spectrum_list(j) = [];
            else
                % Trim buffers
                obj.Peptides(i).spectrum_list(j).peptidoform_list_str(num+1:end) = [];
                obj.Peptides(i).spectrum_list(j).peptidoform_list_abun(num+1:end) = [];
                if isfield(obj.Peptides(i).spectrum_list(j), 'peptidoform_list_support_freq')
                    obj.Peptides(i).spectrum_list(j).peptidoform_list_support_freq(num+1:end) = [];
                end
                if isfield(obj.Peptides(i).spectrum_list(j), 'peptidoform_list_vif')
                    obj.Peptides(i).spectrum_list(j).peptidoform_list_vif(num+1:end) = [];
                end
                if isfield(obj.Peptides(i).spectrum_list(j), 'peptidoform_list_abundance_mad')
                    obj.Peptides(i).spectrum_list(j).peptidoform_list_abundance_mad(num+1:end) = [];
                end
            end
        end
        if isempty(obj.Peptides(i).spectrum_list)
            obj.Peptides(i) = [];
        end
    end
end

% Reset indices just in case
obj.CurrentPeptideIdx = length(obj.Peptides);

obj.rebuildPeptideIndexMap();
end