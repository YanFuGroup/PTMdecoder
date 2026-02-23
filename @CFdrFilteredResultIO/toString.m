function out = toString(value)
% Convert a value to string for table output.
% Input:
%   value (any)
%       Value to convert for output formatting.
% Output:
%   out (1 x 1 char)
%       String representation of the value.
%
% Numeric values are formatted with a compact precision to preserve
% Mascot/FDR table fidelity while avoiding excessive digits.

if isstring(value)
    out = char(value);
elseif ischar(value)
    out = value;
elseif isnumeric(value)
    if isempty(value)
        out = '';
    elseif isscalar(value)
        out = num2str(value, '%.6f');
    else
        out = num2str(value);
    end
else
    out = char(string(value));
end
end
