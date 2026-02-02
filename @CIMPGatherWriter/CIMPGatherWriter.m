classdef CIMPGatherWriter
    % CIMPGatherWriter
    % Writer for IMP gather quantification outputs.
    methods (Static)
        start_new_run(output_path)

        write_imp_group_block(output_path, protein_name_pos, imp_records)

        write_imp_records(fid, imp_records)

        write_header(fid)

        write_protein_start_position_line(fid, prot_names_pos)
    end
end
