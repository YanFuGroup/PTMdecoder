classdef CProteinPeptidePairFileReader
    % Utility for reading protein-peptide pair files.

    methods (Static)
        function [prot_list, peptide_list] = readPairFile(file_path)
            % Read protein-peptide pairs from a tab-delimited file.
            % Input:
            %   file_path (1 x 1 char/string)
            %       Pair file path. Header must be: protein<TAB>peptide.
            % Output:
            %   prot_list (N x 1 cell)
            %       Protein names from first column.
            %   peptide_list (N x 1 cell)
            %       Peptide sequences from second column.

            fid = fopen(file_path, 'r');
            if fid <= 0
                error('CProteinPeptidePairFileReader:OpenPairFileFailed', ...
                    'Failed to open protein-peptide pair file: %s', file_path);
            end

            header_line = fgetl(fid);
            if ~ischar(header_line)
                fclose(fid);
                error('CProteinPeptidePairFileReader:EmptyPairFile', ...
                    'Protein-peptide pair file is empty: %s', file_path);
            end

            header_tokens = regexp(strtrim(header_line), '\t', 'split');
            if numel(header_tokens) ~= 2 || ...
                    ~strcmpi(strtrim(header_tokens{1}), 'protein') || ...
                    ~strcmpi(strtrim(header_tokens{2}), 'peptide')
                fclose(fid);
                error('CProteinPeptidePairFileReader:InvalidPairHeader', ...
                    'The pair file header must be exactly: protein<TAB>peptide.');
            end

            prot_list = {};
            peptide_list = {};
            line_num = 1;
            while ~feof(fid)
                line_num = line_num + 1;
                str_line = fgetl(fid);
                if ~ischar(str_line)
                    continue;
                end

                str_line = strtrim(str_line);
                if isempty(str_line)
                    continue;
                end

                cols = regexp(str_line, '\t', 'split');
                if numel(cols) ~= 2
                    fclose(fid);
                    error('CProteinPeptidePairFileReader:InvalidPairRow', ...
                        'Line %d in %s must have exactly 2 TAB-separated columns.', line_num, file_path);
                end

                prot_name = strtrim(cols{1});
                pep_name = strtrim(cols{2});
                if isempty(prot_name) || isempty(pep_name)
                    fclose(fid);
                    error('CProteinPeptidePairFileReader:EmptyPairField', ...
                        'Line %d in %s contains an empty protein or peptide field.', line_num, file_path);
                end

                prot_list{end+1,1} = prot_name; %#ok<AGROW>
                peptide_list{end+1,1} = pep_name; %#ok<AGROW>
            end
            fclose(fid);

            if isempty(prot_list)
                error('CProteinPeptidePairFileReader:NoPairsFound', ...
                    'No protein-peptide pairs found in file: %s', file_path);
            end
        end
    end
end
