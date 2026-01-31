function runQuant(obj)
% Execute the normalization peptide calculation for all experiments
% Input:
%   obj (CPepNormalization)
%       Normalization processor instance
    
    fprintf('Starting quantifying normalization peptides...\n');
    
    % Process each experiment
    for i_exp = 1:length(obj.experimentNames)
        experimentName = obj.experimentNames{i_exp};

        fprintf('Processing experiment %d/%d: %s\n', i_exp, ...
            length(obj.experimentNames), experimentName);
        
        % Indexing the ms1 and ms2
        specPath = fullfile(obj.spectra_dir, experimentName);
        ms12DatasetIO = CMS12DatasetIO(specPath, obj.ms1_tolerance);
        ms12DatasetIO.SetMap();

        % Open the FDR filtered result
        filtered_res_file_path = fullfile(obj.result_dir, experimentName, obj.input_file_name);
        fin = fopen(filtered_res_file_path, 'r');
        if fin == -1
            error('Cannot open the FDR filtered result file: "%s"!', filtered_res_file_path);
        end

        % Check and create a new output file
        outputPath = fullfile(obj.result_dir, experimentName);
        obj.createOutputFile(outputPath);

        % Initialize quantification objects
        pep_quant = obj.initializeQuantificationObjects(outputPath, ms12DatasetIO);

        % Process the filtered result file
        fprintf('Reading %s...', experimentName);
        pep_quant = obj.readSearchResult(fin, filtered_res_file_path, ms12DatasetIO, pep_quant);
        fprintf('done.\n');

        % Run quantification
        fprintf('Quantifying %s...', experimentName);
        for i_list = 1:length(obj.peptide_list)
            pep_quant{i_list}.runGather();
        end
        fprintf('done.\n');

        % Close files
        fclose(fin);
    end
    
    fprintf('Quantification of normalization peptides completed.\n');
end