function createOutputFile(obj, outputPath)
    % Create and initialize the output file

    fout = fopen(fullfile(outputPath, obj.output_file_name), 'w');
    if fout <= 0
        error('Cannot open the report file %s', ...
            fullfile(outputPath, obj.output_file_name));
    end
    
    fprintf(fout, 'Protein_name,Peptide_start_position_on_protein;\n');
    fprintf(fout, '*\tPeptide\tCharge\tDataset\tMass_center\tLow_mass_bound\tHigh_mass_bound\tPeak_area\n');
    fprintf(fout, '@\tRT_start\tRT_end\tProportion\tCheck_label\n');
    fclose(fout);
end