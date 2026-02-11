classdef CFdrFilteredResultIO
    methods (Static)
        function result = read(filePath)
            % Read FDR filtered result table from file
            fin = fopen(filePath, 'r');
            if fin == -1
                error('Cannot open the FDR filtered result file: "%s"!', filePath);
            end
            cleanup = onCleanup(@() fclose(fin));

            fileInfo = dir(filePath);
            progress_printer = CPrintProgress(fileInfo.bytes);

            header = fgetl(fin);
            if ~ischar(header) && ~isstring(header)
                error('Empty FDR filtered result file: "%s"!', filePath);
            end

            entries = CFdrFilteredResult.emptyEntries();
            lineNo = 1;
            while ~feof(fin)
                strline = fgetl(fin);
                lineNo = lineNo + 1;
                progress_printer = progress_printer.update_show(ftell(fin));

                if isempty(strline)
                    continue;
                end

                segment = regexp(strline, '\t| +', 'split');
                segment = segment(~cellfun('isempty', segment));

                if numel(segment) < 14
                    error('Invalid FDR filtered line %d in "%s": %s', lineNo, filePath, strline);
                end

                entry = struct( ...
                    'Site', segment{1}, ...
                    'DatasetName', segment{2}, ...
                    'Scan', segment{3}, ...
                    'Spectrum', segment{4}, ...
                    'Charge', segment{5}, ...
                    'Calc_neutral_pepmass', segment{6}, ...
                    'precursor_neutral_mass', segment{7}, ...
                    'massdiff', segment{8}, ...
                    'num_match_ions', segment{9}, ...
                    'peptide', segment{10}, ...
                    'protein', segment{11}, ...
                    'modification', segment{12}, ...
                    'modificationlocation', segment{13}, ...
                    'Score', segment{14} ...
                );
                entries(end + 1) = entry; %#ok<AGROW>
            end

            progress_printer.last_update();
            result = CFdrFilteredResult(entries);
        end

        function write(resultOrEntries, filename)
            % Write Mascot result table to file
            if isa(resultOrEntries, 'CFdrFilteredResult')
                entries = resultOrEntries.entries;
            else
                entries = resultOrEntries;
            end

            fout = fopen(filename, 'w');
            if fout == -1
                error('Failed to open the file "%s"!', filename);
            end
            cleanup = onCleanup(@() fclose(fout));

            fprintf(fout, ['Site\tDatasetName\tScan\tSpectrum\tCharge\t' ...
                'Calc_neutral_pepmass\tprecursor_neutral_mass\tmassdiff' ...
                '\tnum_match_ions\tpeptide\tprotein\tmodification\t' ...
                'modificationlocation\tScore\n']);

            for i_entry = 1:length(entries)
                e = entries(i_entry);
                fields = { ...
                    e.Site, e.DatasetName, e.Scan, e.Spectrum, e.Charge, ...
                    e.Calc_neutral_pepmass, e.precursor_neutral_mass, e.massdiff, ...
                    e.num_match_ions, e.peptide, e.protein, e.modification, ...
                    e.modificationlocation, e.Score ...
                };
                fields = cellfun(@CFdrFilteredResultIO.toString, fields, 'UniformOutput', false);
                fprintf(fout, '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n', fields{:});
            end
        end

        function out = toString(value)
            if isstring(value)
                out = char(value);
            elseif ischar(value)
                out = value;
            elseif isnumeric(value)
                if isempty(value)
                    out = '';
                elseif isscalar(value)
                    out = num2str(value, '%.15g');
                else
                    out = num2str(value);
                end
            else
                out = char(string(value));
            end
        end
    end
end
