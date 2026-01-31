% WRITE_REPORT_MSMS_TOP1 Writes a report in the specified format.
%
%   WRITE_REPORT_MSMS_TOP1(filename, result) writes a report to the file
%   specified by filename using the data in result. The report is written
%   in a specific format where each peptide and its associated spectra are
%   listed.
%
%   INPUT:
%       result  - A struct array containing the following fields:
%           .peptide              - A string representing the peptide sequence.
%           .modification         - A string representing the modification.
%           .modificationlocation - A string representing the modification location.
%           .DatasetName          - A string representing the dataset name.
%           .Spectrum             - A string representing the spectrum.
%       filename    - A string specifying the name of the output file.
%
%   The function writes the data to the specified file in the following format:
%       - Each new peptide starts with a line beginning with 'P' followed by the peptide sequence.
%       - Each spectrum associated with the peptide is listed on a new line starting with 'S'.
%       - The modified peptide sequence and a fixed value of 1.000000 are also written.
%
%   Example:
%       result(1).peptide = 'PEPTIDEK';
%       result(1).modification = 'Carbamidomethyl,Oxidation';
%       result(1).modificationlocation = '2,5';
%       result(1).DatasetName = 'Dataset1';
%       result(1).Spectrum = 'Spectrum1';
%       write_report_msms_top1('output.txt', result);
function write_report_msms_top1(result, filename)
% Input:
%   result (1 x N struct)
%       fields: peptide, modification, modificationlocation, DatasetName, Spectrum
%   filename (1 x 1 char/string)
%       output file path

fout = fopen(filename,'w');
if 0 >= fout
    error(['Filed to open the result file report_msms_top1.txt at: "',filename,'"!']);
end

cur_pep = ''; % temp of current peptide, if changed, new pep_spec set is started
for i_out = 1:length(result)
    if ~isequal(result(i_out).peptide,cur_pep)
        if ~isempty(cur_pep)
            fprintf(fout,'\n\n');
        end
        cur_pep = result(i_out).peptide;
        fprintf(fout,['P\t', cur_pep, '\n']);
    end
    mod_pep = get_mod_pep(result(i_out).peptide, result(i_out).modification, result(i_out).modificationlocation);
    fprintf(fout,['S\t', result(i_out).DatasetName,'\t',...
        result(i_out).Spectrum,'\n',...
        mod_pep,'\t1.000000\n']);
end

fclose(fout);
end



function mod_pep = get_mod_pep(pep, mod_name, mod_pos)
% Get the modified peptide string
% Input:
%   pep (1 x 1 char/string): peptide sequence
%   mod_name (1 x 1 char/string): modification name(s)
%   mod_pos (1 x 1 char/string): modification position(s)
% Output:
%   mod_pep (1 x 1 char/string): modified peptide string

mod_pep = ['_', pep, '_'];
if ~isequal(mod_name,'-')
    mod_names = strsplit(mod_name, ',');
    for i = 1:length(mod_names)
        mod_names{i} = strsplit(mod_names{i}, ' ');
        mod_names{i} = mod_names{i}{1};
    end
    mod_positions = strsplit(mod_pos, ',');
    for i = length(mod_names):-1:1
        mod_pep = [mod_pep(1:str2double(mod_positions{i})+1), '{' , ...
            mod_names{i}, '}', mod_pep(str2double(mod_positions{i})+2:end)];
    end
end
end