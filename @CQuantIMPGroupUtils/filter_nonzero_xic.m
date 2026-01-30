function [idxNonZero, auxic, rt_bound, varargout] = filter_nonzero_xic(auxic, rt_bound, varargin)
% Filter outputs by non-zero XIC area.
%
% Inputs:
%   auxic (K x 1 double) area
%       Area under curve for each IMP
%   rt_bound (K x P struct)
%       RT bounds per IMP per peak, fields: .start/.end (minutes)
%   varargin (cell)
%       Additional arrays to filter using idxNonZero (each K x ...)
%
% Outputs:
%   idxNonZero (M x 1 double)
%       Indices of non-zero area rows
%   auxic (M x 1 double)
%       Filtered auxic
%   rt_bound (M x P struct)
%       RT bounds per IMP per peak, fields: .start/.end (minutes)
%   varargout
%       Filtered extra inputs with the same row selection

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
