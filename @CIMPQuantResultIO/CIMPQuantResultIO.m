classdef CIMPQuantResultIO
    % Read/write utilities for IMP quantification reports

    methods (Static)
        obj = read(path)

        write(report, path)

        % Legacy writer helpers (shared by CIMPGatherWriter)
        write_header(fid)

        write_imp_group_block(output_path, protein_name_pos, imp_records)

        write_imp_records(fid, imp_records)

        write_protein_start_position_line(fid, prot_names_pos)
    end

    methods (Static, Access=private)
        protein_name_pos = parse_protein_line(strline)

        rec = parse_record_line(strline)

        peak = parse_rt_peak_line(strline)
    end
end
