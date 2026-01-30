function is_reserved = hasMinRows(~, ratio_matrix, min_rows)
% Check whether the XIC peaks have at least a minimum number of PSMs.
% Inputs:
%   ratio_matrix (N x K double)
%       quantification matrix derived from the MSMS deconvolution
%   min_rows (1 x 1 double/int)
%       minimum number of rows required
% Outputs:
%   is_reserved (1 x 1 logical)
%       whether the matrix has at least min_rows rows

if nargin < 3
    min_rows = 1; % Default minimum number of rows
end

is_reserved = size(ratio_matrix, 1) >= min_rows;

end