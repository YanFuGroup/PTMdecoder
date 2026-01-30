function [idxNonZero, auxic, rt_bound, varargout] = filter_nonzero_xic(auxic, rt_bound, varargin)
% Filter outputs by non-zero XIC area.
%
% Inputs:
%   auxic: area under curve for each IMP
%   rt_bound: retention time bounds (struct array or numeric)
%   varargin: additional arrays to filter using idxNonZero
%
% Outputs:
%   idxNonZero: indices of non-zero area rows
%   auxic: filtered auxic
%   rt_bound: filtered rt_bound
%   varargout: filtered extra inputs

idxNonZero = find(auxic(:,1) ~= 0);
auxic = auxic(idxNonZero,:);
rt_bound = rt_bound(idxNonZero,:);

varargout = cell(size(varargin));
for k = 1:numel(varargin)
    val = varargin{k};
    if isempty(val)
        varargout{k} = val;
    else
        varargout{k} = val(idxNonZero,:);
    end
end
end
