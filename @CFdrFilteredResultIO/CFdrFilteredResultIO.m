classdef CFdrFilteredResultIO
    methods (Static)
        % Read FDR filtered result table from file.
        result = read(filePath)

        % Write a result table to file.
        write(filteredResults, filename)

        % Convert value to string for writing.
        out = toString(value)
    end
end
