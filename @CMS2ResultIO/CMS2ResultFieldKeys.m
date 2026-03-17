function keys = CMS2ResultFieldKeys()
% CMS2ResultFieldKeys - Centralized field keys for report_msms named fields
% Output:
%   keys (1 x 1 struct)
%       keys.spectrum.jaccard : spectrum-level stability key
%       keys.spectrum.condX   : reserved spectrum-level key
%       keys.peptidoform.support : peptidoform-level support key
%       keys.peptidoform.mad : peptidoform-level abundance MAD key

keys = struct();
keys.spectrum = struct();
keys.peptidoform = struct();

keys.spectrum.jaccard = 'jaccard';
keys.spectrum.condX = 'condX';
keys.peptidoform.support = 'support';
keys.peptidoform.mad = 'mad';
end
