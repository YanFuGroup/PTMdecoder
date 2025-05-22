function write_peptide_spectra_list_file(result, filename)

fout = fopen(filename,'w');
if 0 >= fout
    error(['Filed to open dataset file "',filename,'"!']);
end

cur_pep=''; % temp of current peptide, if changed, new pep_spec set is started
for i_out = 1:length(result)
    if ~isequal(result(i_out).peptide,cur_pep)
        cur_pep = result(i_out).peptide;
        fprintf(fout,[cur_pep,'\n']);
    end
    fprintf(fout,[result(i_out).DatasetName,'\t',...
        result(i_out).Spectrum,'\n']);
end

% close those files
fclose(fout);
end