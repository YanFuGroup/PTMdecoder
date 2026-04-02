function obj = site_level_dataset_summary(obj)
% Prepare site-by-dataset matrix data from peptide-level result.
% Inputs:
%   obj (CSiteLevelDatasetSummary)
%       dataset-level site summarizer instance
% Outputs:
%   obj (CSiteLevelDatasetSummary)
%       updated with full site-by-dataset aggregation

fin = fopen(obj.m_input_path);
if fin < 0
    error(['Can not open the peptide level result file: "', obj.m_input_path, '"']);
end

CLogger.info(['[CSiteLevelDatasetSummary:site_level_dataset_summary] Start summarizing. ', ...
    'input=%s'], obj.m_input_path);

% -------------------------------------------------------------------------
% The peptide-level result file contains three fixed-format header lines,
% then a mixed stream of protein-context lines and peptide-detail lines.
% We build a two-level map:
%   site_dataset_sum(site_name) -> dataset_sum_map(dataset_name) = AUC sum
% and track all observed datasets for final matrix column ordering.
% -------------------------------------------------------------------------

% Keep compatibility with peptide-level file format:
% line 1: protein header, line 2: peptide header, line 3: xic peak header
fgetl(fin);
fgetl(fin);
fgetl(fin);

site_dataset_sum = containers.Map('KeyType', 'char', 'ValueType', 'any');
dataset_seen = containers.Map('KeyType', 'char', 'ValueType', 'logical');

line_count = 0;
total_protein_lines = 0;
matched_protein_lines = 0;
unmatched_protein_lines = 0;
peptide_line_count = 0;
malformed_peptide_line_count = 0;
invalid_auc_line_count = 0;
uninterested_peptide_lines = 0;
invalid_protein_context_segment_count = 0;

max_column_idx = max([obj.m_column_idxs.icol_seq, obj.m_column_idxs.icol_dataset, obj.m_column_idxs.icol_auc]);

selected_protein_contexts = struct('protein_name', {}, 'start_pos', {});

% -------------------------------------------------------------------------
% Streaming parser over the whole file.
% - Peptide line: starts with '*', carries sequence / dataset / AUC.
% - Protein line: non-'*' and non-'@', refreshes current protein context.
% The current protein context is then used by following peptide lines until
% the next protein line appears.
% -------------------------------------------------------------------------
while ~feof(fin)
    strline = fgetl(fin);
    if ~ischar(strline) || isempty(strline)
        continue;
    end

    line_count = line_count + 1;
    if strline(1) == '*'
        peptide_line_count = peptide_line_count + 1;

        % Parse and validate peptide-level columns first.
        segments = strsplit(strline, sprintf('\t'));
        if numel(segments) < max_column_idx
            malformed_peptide_line_count = malformed_peptide_line_count + 1;
            continue;
        end

        dataset_name = strtrim(segments{obj.m_column_idxs.icol_dataset});
        if isempty(dataset_name)
            malformed_peptide_line_count = malformed_peptide_line_count + 1;
            continue;
        end

        if isempty(selected_protein_contexts)
            uninterested_peptide_lines = uninterested_peptide_lines + 1;
            continue;
        end

        auc_value = str2double(strtrim(segments{obj.m_column_idxs.icol_auc}));
        if isnan(auc_value) || ~isfinite(auc_value)
            invalid_auc_line_count = invalid_auc_line_count + 1;
            continue;
        end

        dataset_seen(dataset_name) = true;

        modified_peptides = segments{obj.m_column_idxs.icol_seq};

        % Remove configured tokens that should not affect mod parsing.
        for idx_ig = 1:length(obj.m_ignore_strings)
            modified_peptides = strrep(modified_peptides, obj.m_ignore_strings{idx_ig}, '');
        end

        % Extract all modification blocks in the form {mod_name}.
        mod_str_matches = regexp(modified_peptides, '{(.*?)}', 'tokens');
        if isempty(mod_str_matches)
            continue;
        end

        positions_seq = zeros(1, numel(mod_str_matches));
        positions_str = zeros(1, numel(mod_str_matches));
        start_pos = 0;

        % Iterate each modification occurrence and map it to a site key.
        for i_mod = 1:numel(mod_str_matches)
            mod_name = mod_str_matches{i_mod}{1};
            found_index = strfind(modified_peptides(start_pos + 1:end), ['{' mod_name '}']);
            if isempty(found_index)
                continue;
            end

            if i_mod == 1
                positions_seq(i_mod) = found_index(1) - 2;
            else
                positions_seq(i_mod) = positions_seq(i_mod - 1) + found_index(1) - 1;
            end

            positions_str(i_mod) = start_pos + found_index(1);
            start_pos = positions_str(i_mod) + numel(mod_name) + 1;

            if ~isKey(obj.m_mod_name_abbr, mod_name)
                continue;
            end

            abbr_mod = obj.m_mod_name_abbr(mod_name);
            mod_specificity = modified_peptides(positions_str(i_mod) - 1);

            % Build one protein-prefix string that preserves all proteins from
            % the current protein line, each with its computed site position.
            protein_site_tokens = cell(1, numel(selected_protein_contexts));
            for idx_ctx = 1:numel(selected_protein_contexts)
                mod_prot_pos = selected_protein_contexts(idx_ctx).start_pos + positions_seq(i_mod) - 1;
                if obj.m_site_position_count_initial_m
                    site_pos_curr = mod_prot_pos;
                else
                    site_pos_curr = mod_prot_pos - 1;
                end

                protein_site_tokens{idx_ctx} = [selected_protein_contexts(idx_ctx).protein_name, ...
                    ',', num2str(site_pos_curr), ';'];
            end
            protein_site_prefix = strjoin(protein_site_tokens, '');

            % Keep terminal/residue semantics from the legacy output format.
            if mod_specificity == '_'
                if positions_seq(i_mod) == 0
                    mod_specificity = 'N-term';
                else
                    mod_specificity = 'C-term';
                end
            end
            % Main site name format:
            % [protein_name,pos;protein_name,pos;] [specificity]_[abbr]
            site_name = [protein_site_prefix, ' ', mod_specificity, '_', abbr_mod];

            % Aggregate AUC by (site, dataset).
            if isKey(site_dataset_sum, site_name)
                dataset_sum_map = site_dataset_sum(site_name);
            else
                dataset_sum_map = containers.Map('KeyType', 'char', 'ValueType', 'double');
            end

            if isKey(dataset_sum_map, dataset_name)
                dataset_sum_map(dataset_name) = dataset_sum_map(dataset_name) + auc_value;
            else
                dataset_sum_map(dataset_name) = auc_value;
            end

            site_dataset_sum(site_name) = dataset_sum_map;
        end
    elseif strline(1) ~= '@'
        % Protein-context line: cache all protein/start-position pairs from
        % this line and reuse them for following peptide-detail lines.
        total_protein_lines = total_protein_lines + 1;
        selected_protein_contexts = struct('protein_name', {}, 'start_pos', {});

        segments = strsplit(strline, ';');
        for i_seg = 1:(length(segments) - 1)
            key_value = strsplit(segments{i_seg}, ',');
            if numel(key_value) < 2
                invalid_protein_context_segment_count = invalid_protein_context_segment_count + 1;
                continue;
            end
            protein_name = strtrim(key_value{1});
            start_pos = str2double(strtrim(key_value{2}));
            if isempty(protein_name) || isnan(start_pos) || ~isfinite(start_pos)
                invalid_protein_context_segment_count = invalid_protein_context_segment_count + 1;
                continue;
            end

            new_idx = numel(selected_protein_contexts) + 1;
            selected_protein_contexts(new_idx).protein_name = protein_name;
            selected_protein_contexts(new_idx).start_pos = start_pos;
        end

        if isempty(selected_protein_contexts)
            unmatched_protein_lines = unmatched_protein_lines + 1;
        else
            matched_protein_lines = matched_protein_lines + 1;
        end
    end
end

fclose(fin);

% -------------------------------------------------------------------------
% Finalize outputs: sort dataset/site axes for deterministic downstream use,
% then emit quality diagnostics for skipped or unmatched records.
% -------------------------------------------------------------------------
dataset_names = keys(dataset_seen);
if isempty(dataset_names)
    dataset_names = {};
else
    dataset_names = sort(dataset_names);
end

site_names = keys(site_dataset_sum);
if isempty(site_names)
    site_names = {};
else
    site_names = sort(site_names);
end

obj.m_site_dataset_sum = site_dataset_sum;
obj.m_dataset_names = dataset_names;
obj.m_site_names = site_names;

if malformed_peptide_line_count > 0
    CLogger.warn(['[CSiteLevelDatasetSummary:site_level_dataset_summary] ', ...
        'Malformed peptide lines skipped: %d.'], malformed_peptide_line_count);
end

if invalid_auc_line_count > 0
    CLogger.warn(['[CSiteLevelDatasetSummary:site_level_dataset_summary] ', ...
        'Peptide lines with invalid AUC skipped: %d.'], invalid_auc_line_count);
end

if unmatched_protein_lines > 0
    CLogger.warn(['[CSiteLevelDatasetSummary:site_level_dataset_summary] ', ...
        'Unmatched protein lines found: %d (mapped=%d, total=%d).'], ...
        unmatched_protein_lines, matched_protein_lines, total_protein_lines);
end

if invalid_protein_context_segment_count > 0
    CLogger.warn(['[CSiteLevelDatasetSummary:site_level_dataset_summary] ', ...
        'Invalid protein-context entries skipped: %d.'], invalid_protein_context_segment_count);
end

CLogger.info(['[CSiteLevelDatasetSummary:site_level_dataset_summary] Done. ', ...
    'lines=%d, protein_lines=%d, peptide_lines=%d, datasets=%d, sites=%d, uninterested_peptides=%d.'], ...
    line_count, total_protein_lines, peptide_line_count, numel(dataset_names), numel(site_names), uninterested_peptide_lines);

end
