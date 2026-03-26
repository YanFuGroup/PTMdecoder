function keys = CMS2ResultFieldKeys()
% CMS2ResultFieldKeys - Centralized field keys for report_msms named fields
% Output:
%   keys (1 x 1 struct)
%       keys.spectrum.jaccard : spectrum-level stability key
%       keys.spectrum.vifAll  : spectrum-level max VIF on all IMP columns
%       keys.spectrum.vifReported : spectrum-level max VIF on reported IMP columns
%       keys.peptidoform.support : peptidoform-level support key
%       keys.peptidoform.vif : peptidoform-level VIF key
%       keys.peptidoform.mad : peptidoform-level abundance MAD key

keys = struct();
keys.spectrum = struct();
keys.peptidoform = struct();

keys.spectrum.jaccard = 'jaccard';
keys.spectrum.vifAll = 'vif_all';
keys.spectrum.vifReported = 'vif_reported';
keys.peptidoform.support = 'support';
keys.peptidoform.vif = 'vif';
keys.peptidoform.mad = 'mad';
end
