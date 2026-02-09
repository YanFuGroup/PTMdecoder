classdef CIMPQuantResultIO
    % Read/write utilities for IMP quantification reports

    methods (Static)
        obj = read(path)

        write(report, path)
    end

    methods (Static, Access=private)
        protein_name_pos = parse_protein_line(strline)

        rec = parse_record_line(strline)

        peak = parse_rt_peak_line(strline)
    end
end
