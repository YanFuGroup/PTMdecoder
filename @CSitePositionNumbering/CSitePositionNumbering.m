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
            % regexp otherwise tolerates malformed syntax as literal text.
            previous_state = warning;
            [previous_message, previous_id] = lastwarn;
            cleanup = onCleanup(@() CSitePositionNumbering.restoreWarnings( ...
                previous_state, previous_message, previous_id)); %#ok<NASGU>
            warning('on', 'all');
            warning('off', 'backtrace');
            lastwarn('');
            try
                regexp('>PTMdecoder_regex_validation', expression, 'tokens', 'once', 'warnings');
            catch cause
                error('CSitePositionNumbering:InvalidRegularExpression', ...
                    'Invalid FASTA protein-name regular expression: %s', cause.message);
            end
            [message, id] = lastwarn;
            if startsWith(id, 'MATLAB:REGEXP:')
                error('CSitePositionNumbering:InvalidRegularExpression', ...
                    'Invalid FASTA protein-name regular expression: %s', message);
            end
        end

        function restoreWarnings(state, message, id)
            warning(state);
            lastwarn(message, id);
        end
    end
end
