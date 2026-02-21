function addPeptidoform(obj, peptidoform_str, relative_abundance)
% Input:
%   peptidoform_str (1 x 1 char/string)
%       peptidoform string
%   relative_abundance (1 x 1 double)
%       relative abundance

if obj.CurrentPeptideIdx == 0 || obj.CurrentSpectrumIdx == 0
    error('CMS2Result:NoSpectrum', 'Cannot add peptidoform without peptide and spectrum context.');
end

% Quick access references
currentNum = obj.Peptides(obj.CurrentPeptideIdx).spectrum_list(obj.CurrentSpectrumIdx).peptidoform_num;
newNum = currentNum + 1;

% Update count
obj.Peptides(obj.CurrentPeptideIdx).spectrum_list(obj.CurrentSpectrumIdx).peptidoform_num = newNum;

% Buffer check logic from original code
currentLen = length(obj.Peptides(obj.CurrentPeptideIdx).spectrum_list(obj.CurrentSpectrumIdx).peptidoform_list_abun);
if newNum > currentLen
    % Extend by 50
    % TODO: make the buffer size a parameter?
    obj.Peptides(obj.CurrentPeptideIdx).spectrum_list(obj.CurrentSpectrumIdx).peptidoform_list_str{newNum + 50} = '';
    obj.Peptides(obj.CurrentPeptideIdx).spectrum_list(obj.CurrentSpectrumIdx).peptidoform_list_abun(newNum + 50) = 0;
end

obj.Peptides(obj.CurrentPeptideIdx).spectrum_list(obj.CurrentSpectrumIdx).peptidoform_list_str{newNum} = peptidoform_str;
obj.Peptides(obj.CurrentPeptideIdx).spectrum_list(obj.CurrentSpectrumIdx).peptidoform_list_abun(newNum) = relative_abundance;
end