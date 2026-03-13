function jaccard = computeJaccardIndex(maskA, maskB)
% computeJaccardIndex - Compute Jaccard index between two logical masks
% Inputs:
%   maskA (N x 1 logical or 1 x N logical)
%       First logical membership mask.
%   maskB (N x 1 logical or 1 x N logical)
%       Second logical membership mask.
% Outputs:
%   jaccard (1 x 1 double)
%       Jaccard index. If both masks are empty sets, return 0.
%       Attention: Jaccard index is not defined for empty sets, but we return 0 to avoid NaN and to reflect no shared membership.

if ~islogical(maskA) || ~islogical(maskB)
    CLogger.error('[CMS2QuantSolver:InvalidMaskType] maskA and maskB must be logical arrays.');
end

maskA = maskA(:);
maskB = maskB(:);

if numel(maskA) ~= numel(maskB)
    CLogger.error(['[CMS2QuantSolver:InvalidMaskSize] ', ...
        'maskA and maskB must have the same number of elements.']);
end

intersection_count = sum(maskA & maskB);
union_count = sum(maskA | maskB);

if union_count == 0
    jaccard = 0;
    return;
end

jaccard = intersection_count / union_count;
end
