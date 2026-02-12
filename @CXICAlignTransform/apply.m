function rt_out = apply(model, rt_in)
% Apply a fitted RT alignment transform.
% Input:
%   model (struct)
%       Transform model from CXICAlignTransform.fit
%   rt_in (double array)
%       Input RTs (minutes)
% Output:
%   rt_out (double array)
%       Aligned RTs (minutes)

rt_out = model.slope .* rt_in + model.intercept;
if ~isempty(model.bin_centers)
    % Avoid excessive extrapolation by clamping to bin range
    rt_query = rt_in;
    if ~isempty(model.bin_min) && ~isempty(model.bin_max)
        rt_query = min(max(rt_query, model.bin_min), model.bin_max);
    end
    rt_out = rt_out + interp1(model.bin_centers, model.bin_offsets, rt_query, 'linear', 'extrap');
end
end
