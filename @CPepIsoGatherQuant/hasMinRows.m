function is_reserved = hasMinRows(obj, ratio_matrix, min_rows)
% Check whether the XIC peaks with at least THREE PSMs.
% Inputs:
%   ratio_matrix: the quantification matrix derived from the MSMS deconvolution
% Outputs:
%   is_reserved: a logical value indicating whether each XIC peak is reserved

if nargin < 2
    min_rows = 1; % Default minimum number of rows
end

is_reserved = size(ratio_matrix, 1) >= min_rows;

end