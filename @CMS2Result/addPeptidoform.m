function addPeptidoform(obj, peptidoform_str, relative_abundance, varargin)
% Input:
%   peptidoform_str (1 x 1 char/string)
%       peptidoform string
%   relative_abundance (1 x 1 double)
%       relative abundance
%   varargin{1} (optional, 1 x 1 struct)
%       optional support frequency metadata
%       - support_frequency (1 x 1 double)
%       - abundance_mad (1 x 1 double, reserved for future use)

if obj.CurrentPeptideIdx == 0 || obj.CurrentSpectrumIdx == 0
    error('CMS2Result:NoSpectrum', 'Cannot add peptidoform without peptide and spectrum context.');
end

support_frequency = NaN;
if nargin >= 4 && ~isempty(varargin{1})
    supportMeta = varargin{1};
    if ~isstruct(supportMeta)
        error('CMS2Result:InvalidSupportFrequency', ...
            'Optional support metadata must be a struct.');
    end
    if isfield(supportMeta, 'support_frequency')
        support_frequency = supportMeta.support_frequency;
    end
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
    obj.Peptides(obj.CurrentPeptideIdx).spectrum_list(obj.CurrentSpectrumIdx).peptidoform_list_support_freq(newNum + 50) = NaN;
end

obj.Peptides(obj.CurrentPeptideIdx).spectrum_list(obj.CurrentSpectrumIdx).peptidoform_list_str{newNum} = peptidoform_str;
obj.Peptides(obj.CurrentPeptideIdx).spectrum_list(obj.CurrentSpectrumIdx).peptidoform_list_abun(newNum) = relative_abundance;
obj.Peptides(obj.CurrentPeptideIdx).spectrum_list(obj.CurrentSpectrumIdx).peptidoform_list_support_freq(newNum) = support_frequency;
end