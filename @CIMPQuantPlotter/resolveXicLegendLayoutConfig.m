function layout = resolveXicLegendLayoutConfig(xic_layout)
% Resolve an optional XIC layout override against the default layout.
% Input:
%   xic_layout (struct or [], optional)
%       partial layout override supplied by code API callers
% Output:
%   layout (struct)
%       complete layout with derived horizontal geometry recalculated

if nargin < 1 || isempty(xic_layout)
    layout = CIMPQuantPlotter.getXicLegendLayoutConfig();
    return;
end

if ~isstruct(xic_layout) || ~isscalar(xic_layout)
    error('CIMPQuantPlotter:InvalidXicLayoutOverride', ...
        'xic_layout must be a scalar struct or empty.');
end

layout = CIMPQuantPlotter.getXicLegendLayoutConfig();
override_fields = fieldnames(xic_layout);
validateNoDerivedFieldOverrides(override_fields);
for idx_field = 1:numel(override_fields)
    field_name = override_fields{idx_field};
    layout.(field_name) = xic_layout.(field_name);
end

layout = recompute_horizontal_geometry(layout);
end

function layout = recompute_horizontal_geometry(layout)
% Keep axes_width_px and legend_max_width_px coherent after overrides.
available_width_px = layout.figure_width_px ...
    - layout.left_margin_px - layout.axes_legend_gap_px ...
    - layout.right_margin_px - layout.legend_print_safety_px;

if isfield(layout, 'axes_width_fraction') && ~isempty(layout.axes_width_fraction)
    layout.axes_width_px = round(available_width_px * layout.axes_width_fraction);
    layout.legend_max_width_px = available_width_px - layout.axes_width_px;
    return;
end

layout.legend_max_width_px = available_width_px - layout.axes_width_px;
end

function validateNoDerivedFieldOverrides(override_fields)
% Reject output-only geometry fields that must stay derived from source layout inputs.
derived_fields = {'axes_width_px', 'legend_max_width_px'};
invalid_mask = ismember(derived_fields, override_fields);
if any(invalid_mask)
    error('CIMPQuantPlotter:InvalidXicLayoutOverride', ...
        ['xic_layout must not override derived fields: %s. ', ...
        'Use figure_width_px and axes_width_fraction instead.'], ...
        strjoin(derived_fields(invalid_mask), ', '));
end
end
