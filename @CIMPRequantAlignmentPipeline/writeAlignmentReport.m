function writeAlignmentReport(~, report, output_path)
% Write alignment report to a text file.
% Input:
%   report (struct)
%       Alignment report with pair stats and summary
%   output_path (1 x 1 char/string)
%       Output path for the report
% Output:
%   (none)

fout = fopen(output_path, 'w');
if fout == -1
    error('Failed to open alignment report file: %s', output_path);
end
cleanup = onCleanup(@() fclose(fout));

fprintf(fout, 'XIC alignment report\n');
fprintf(fout, 'Total anchors: %d\n', report.total_anchors);

if isfield(report, 'pairs') && ~isempty(report.pairs)
    fprintf(fout, '\nPair models:\n');
    fprintf(fout, 'ref_raw\ttarget_raw\tanchors\tslope\tintercept\tmedian_abs_resid\n');
    for idx = 1:numel(report.pairs)
        p = report.pairs(idx);
        fprintf(fout, '%s\t%s\t%d\t%.6f\t%.6f\t%.6f\n', ...
            p.ref_raw, p.target_raw, p.num_anchors, p.slope, p.intercept, p.median_abs_resid);
    end
end

if isfield(report, 'summary')
    s = report.summary;
    fprintf(fout, '\nSummary:\n');
    fprintf(fout, 'aligned_groups\t%d\n', s.num_groups);
    fprintf(fout, 'aligned_imps\t%d\n', s.num_aligned);
    if isfield(s, 'num_missing_pair')
        fprintf(fout, 'missing_pairs\t%d\n', s.num_missing_pair);
    elseif isfield(s, 'num_missing_ref')
        fprintf(fout, 'missing_ref_imps\t%d\n', s.num_missing_ref);
    end
end
end
