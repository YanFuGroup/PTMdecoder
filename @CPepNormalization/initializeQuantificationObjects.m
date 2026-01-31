function pep_quant = initializeQuantificationObjects(obj, outputPath, ms12DatasetIO)
    % Initialize quantification objects for each peptide
    % Input:
    %   obj (CPepNormalization)
    %       Normalization processor instance
    %   outputPath (1 x 1 char/string)
    %       Output directory path
    %   ms12DatasetIO (object)
    %       MS1/MS2 dataset IO instance
    % Output:
    %   pep_quant (K x 1 cell)
    %       Quantification objects per target peptide
    
    pep_quant = cell(length(obj.peptide_list), 1);
    for i_list = 1:length(obj.peptide_list)
        pep_quant{i_list} = CPepIsoGatherQuant({obj.prot_list{i_list}, -1}, ...
            ms12DatasetIO, obj.resFilterThres, obj.ms1_tolerance, obj.alpha, ...
            fullfile(outputPath, obj.output_file_name));
    end
end