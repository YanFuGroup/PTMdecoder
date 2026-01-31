function [imp_idx_nonzero, area_imp_final, rt_bound, varargout] = filter_nonzero_xic(area_imp_final, rt_bound, varargin)
% Filter outputs by non-zero XIC area.
%
% Inputs:
%   area_imp_final (K x 1 double) area
%       Area under curve for each IMP
%   rt_bound (K x P struct)
%       RT bounds per IMP per peak, fields: .start/.end (minutes)
%   varargin (cell)
%       Additional arrays to filter using imp_idx_nonzero (each K x ...)
%
% Outputs:
%   imp_idx_nonzero (M x 1 double)
%       Indices of non-zero area rows
%   area_imp_final (M x 1 double)
%       Filtered area_imp_final
%   rt_bound (M x P struct)
%       RT bounds per IMP per peak, fields: .start/.end (minutes)
%   varargout
%       Filtered extra inputs with the same row selection

imp_idx_nonzero = find(area_imp_final(:,1) ~= 0);
area_imp_final = area_imp_final(imp_idx_nonzero,:);
rt_bound = rt_bound(imp_idx_nonzero,:);

varargout = cell(size(varargin));
for k = 1:numel(varargin)
    val = varargin{k};
    if isempty(val)
        varargout{k} = val;
    else
        varargout{k} = val(imp_idx_nonzero,:);
    end
end
end
