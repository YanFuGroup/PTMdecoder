classdef CSitePositionNumbering < handle
    % Shared, sequence-backed coordinate offset for site report writers.
    % This changes reporting coordinates only; it does not infer M cleavage.
    %
    % Configuration:
    %   count_initial_m is supplied by site_position_count_initial_m.
    %   true returns zero offsets without reading or validating FASTA settings;
    %   this is the matrix and traceability default. false subtracts one only
    %   for a full protein sequence starting with M (case-insensitive), and is
    %   the site-list default. Supply fasta_file_path and regular_express from
    %   the upstream protein mapping; lookup uses the original protein ID.
    %   The expression is passed unchanged to the FASTA reader. Validation uses
    %   regexp syntax, allowing literal pipes in character classes but rejecting
    %   malformed expressions.
    %
    % Required sequence evidence when false:
    %   Omitted FASTA is a configuration error. Referenced IDs must exist and
    %   have nonempty sequences; missing evidence never falls back to zero.
    %   Unavailable/unreadable files, missing or invalid name rules, and
    %   conflicting sequences for one ID also raise errors.
    %   Identical duplicates are accepted. Strict duplicate checking is opt-in
    %   for these reports; other CFastaReader callers retain their defaults.
    %
    % Report integration:
    %   CSiteLevelDatasetSummary and CTraceabilityReportService query each
    %   protein independently. CSiteLevelSummary queries its first mapped
    %   protein and uses unnumbered terminal labels. Matrix/trace N- and C-term
    %   coordinates lie before and after the peptide, respectively. Excluding
    %   initial M can yield M0 or N-term -1 for a peptide starting at protein
    %   position 1; this class does not clamp or reinterpret those coordinates.

    properties (Access = private)
        count_initial_m
        protein_sequences
    end

    methods
        function obj = CSitePositionNumbering(count_initial_m, fasta_file_path, regular_express)
            % Load the same FASTA/name rule used to produce the peptide report.
            % Counting M needs no sequence evidence and performs no FASTA IO.
            obj.count_initial_m = count_initial_m;
            obj.protein_sequences = containers.Map('KeyType', 'char', 'ValueType', 'char');
            if count_initial_m
                return;
            end
            if isempty(fasta_file_path) || strlength(strtrim(string(fasta_file_path))) == 0
                error('CSitePositionNumbering:MissingFasta', ...
                    'fasta_file_path is required when site_position_count_initial_m=false.');
            end
            if ~isfile(fasta_file_path)
                error('CSitePositionNumbering:UnavailableFasta', ...
                    'Configured FASTA file is unavailable: %s', fasta_file_path);
            end
            if strlength(string(regular_express)) == 0
                error('CSitePositionNumbering:InvalidRegularExpression', ...
                    'regular_express is required when fasta_file_path is provided for site numbering.');
            end
            CSitePositionNumbering.validateRegularExpression(regular_express);
            reader = CFastaReader(fasta_file_path, regular_express);
            % Strict duplicate checking is opt-in; search callers are unchanged.
            obj.protein_sequences = reader.read(true);
        end

        function offset = getOffset(obj, protein_name)
            % Subtract one only for a known initial M, using the original ID.
            % Unknown and empty sequences are errors when excluding initial M.
            offset = 0;
            if obj.count_initial_m
                return;
            end
            if ~isKey(obj.protein_sequences, protein_name)
                error('CSitePositionNumbering:MissingProteinSequence', ...
                    ['Protein "%s" is absent from the configured FASTA; ', ...
                    'a complete protein sequence is required when site_position_count_initial_m=false.'], ...
                    protein_name);
            end
            sequence = obj.protein_sequences(protein_name);
            if isempty(sequence)
                error('CSitePositionNumbering:EmptyProteinSequence', ...
                    ['Protein "%s" has an empty sequence in the configured FASTA; ', ...
                    'a complete protein sequence is required when site_position_count_initial_m=false.'], ...
                    protein_name);
            end
            offset = double(strcmpi(sequence(1), 'M'));
        end
    end

    methods (Static, Access = private)
        function validateRegularExpression(expression)
            % Use the reader's regex dialect, including backreferences.
            previous_state = warning;
            previous_backtrace = warning('query', 'backtrace');
            previous_verbose = warning('query', 'verbose');
            [previous_message, previous_id] = lastwarn;
            cleanup = onCleanup(@() CSitePositionNumbering.restoreWarnings( ...
                previous_state, previous_message, previous_id, ...
                previous_backtrace, previous_verbose)); %#ok<NASGU>
            warning('on', 'all');
            warning('off', 'backtrace');
            warning('on', 'verbose');
            lastwarn('');
            try
                diagnostics = evalc(['regexp(''>PTMdecoder_regex_validation'', ', ...
                    'expression, ''tokens'', ''once'', ''warnings'');']);
            catch cause
                error('CSitePositionNumbering:InvalidRegularExpression', ...
                    'Invalid FASTA protein-name regular expression: %s', cause.message);
            end
            % Inspect every diagnostic: lastwarn can hide earlier syntax errors.
            % Remove only MATLAB's warning-off hyperlinks, preserving expression
            % text such as named groups. Verbose IDs do not depend on locale.
            diagnostics = regexprep(diagnostics, ...
                '<a href="matlab:warning off [^"]*">|</a>', '');
            diagnostic_ids = regexp(diagnostics, ...
                'warning off (MATLAB:REGEXP:[A-Za-z0-9_]+)', 'tokens');
            [~, last_id] = lastwarn;
            if isempty(diagnostic_ids) && ~startsWith(last_id, 'MATLAB:REGEXP:')
                return;
            end
            if ~isempty(diagnostic_ids)
                ids = cellfun(@(x) x{1}, diagnostic_ids, 'UniformOutput', false);
                if all(strcmp(ids, 'MATLAB:REGEXP:ErrorLiteral'))
                    % ErrorLiteral also describes malformed ')' or '**'. Only
                    % exempt diagnostics pointing to a pipe inside a character
                    % class; do not rewrite the expression or stop parsing early.
                    expression_text = char(expression);
                    markers = regexp(diagnostics, ['\r?\n\t', ...
                        regexptranslate('escape', expression_text), ...
                        '\r?\n\t( *)\S'], 'tokens');
                    positions = cellfun(@(x) numel(x{1}) + 1, markers);
                    allowed = CSitePositionNumbering.literalPipePositions(expression_text);
                    if numel(positions) == numel(ids) && ...
                            all(positions <= numel(allowed)) && all(allowed(positions))
                        return;
                    end
                end
                last_id = ids{1};
            end
            error('CSitePositionNumbering:InvalidRegularExpression', ...
                'Invalid FASTA protein-name regular expression "%s" (%s).', ...
                expression, last_id);
        end

        function positions = literalPipePositions(expression)
            % Classify diagnostic locations, leaving syntax validation to regexp.
            positions = false(size(expression));
            in_class = false;
            first_member = false;
            idx = 1;
            while idx <= numel(expression)
                ch = expression(idx);
                if ch == '\'
                    idx = idx + 2;
                    first_member = false;
                    continue;
                end
                if ~in_class && ch == '['
                    in_class = true;
                    first_member = true;
                elseif in_class
                    if ch == ']' && ~first_member
                        in_class = false;
                    elseif ch == '|'
                        positions(idx) = true;
                    end
                    if ~(first_member && ch == '^' && expression(idx - 1) == '[')
                        first_member = false;
                    end
                end
                idx = idx + 1;
            end
        end

        function restoreWarnings(state, message, id, backtrace, verbose)
            warning(state);
            warning(backtrace);
            warning(verbose);
            lastwarn(message, id);
        end
    end
end
