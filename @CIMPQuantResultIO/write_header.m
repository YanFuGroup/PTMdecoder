function write_header(fid)
% write_header
% Write file header for IMP quantification output.
% input:
%   fid (1 x 1 double)
%       file identifier for the output file

fprintf(fid,'Protein_name,Peptide_start_position_on_protein;\n');
fprintf(fid,'*\tIMP\tCharge\tDataset\tMass_center\tLow_mass_bound\tHigh_mass_bound\tPeak_area\n');
fprintf(fid,'@\tRT_start\tRT_end\tProportion\tCheck_label\n');
end
