classdef CSitePositionNumbering < handle
    % Shared, sequence-backed coordinate offset for site report writers.
    % This changes reporting coordinates only; it does not infer M cleavage.

    properties (Access = private)
        count_initial_m
        protein_sequences
        missing_sequence_warnings
    end

    methods
        function obj = CSitePositionNumbering(count_initial_m, fasta_file_path, regular_express)
            % Load the same FASTA/name rule used to produce the peptide report.
            % Counting M needs no sequence evidence and performs no FASTA IO.
            obj.count_initial_m = count_initial_m;
            obj.protein_sequences = containers.Map('KeyType', 'char', 'ValueType', 'char');
            obj.missing_sequence_warnings = containers.Map('KeyType', 'char', 'ValueType', 'logical');
            if count_initial_m || strlength(string(fasta_file_path)) == 0
                return;
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
            % Unknown and empty sequences retain coordinates with one warning
            % per protein per numbering instance, even across multiple sites.
            offset = 0;
            if obj.count_initial_m
                return;
            end
            if isKey(obj.protein_sequences, protein_name)
                sequence = obj.protein_sequences(protein_name);
                if ~isempty(sequence)
                    offset = double(strcmpi(sequence(1), 'M'));
                    return;
                end
            end
            if ~isKey(obj.missing_sequence_warnings, protein_name)
                CLogger.warn(['[CSitePositionNumbering:MissingProteinSequence] ', ...
                    'No sequence evidence for protein "%s"; retaining original coordinates ', ...
                    'with site_position_count_initial_m=false.'], protein_name);
                obj.missing_sequence_warnings(protein_name) = true;
            end
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
