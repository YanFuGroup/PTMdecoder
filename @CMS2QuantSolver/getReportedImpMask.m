function reported_imp_mask = getReportedImpMask(abundance, relative_threshold)
% getReportedImpMask - Build reported IMP mask from abundance values
% Inputs:
%   abundance (N x 1 double or 1 x N double)
%       Relative abundance vector.
%   relative_threshold (1 x 1 double)
%       Relative threshold factor.
% Outputs:
%   reported_imp_mask (N x 1 logical)
%       True for IMP entries that satisfy abundance >= tau.

if nargin < 2
	CLogger.error('[CMS2QuantSolver:MissingRelativeThreshold] relative_threshold is required.');
end

if ~isnumeric(abundance)
	CLogger.error('[CMS2QuantSolver:InvalidAbundance] abundance must be a numeric vector.');
end
if ~isscalar(relative_threshold) || ~isnumeric(relative_threshold) || relative_threshold < 0
	CLogger.error(['[CMS2QuantSolver:InvalidRelativeThreshold] ', ...
		'relative_threshold must be a non-negative numeric scalar.']);
end

abundance = abundance(:);

if isempty(abundance)
	reported_imp_mask = false(0, 1);
	return;
end

max_abundance = max(abundance);

if max_abundance <= 0
	reported_imp_mask = false(size(abundance));
	return;
end

reported_imp_mask = abundance >= relative_threshold * max_abundance;
end
