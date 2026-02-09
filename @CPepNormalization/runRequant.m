function runRequant(obj, checked_pep_file_path)
% Execute the normalization peptide calculation for all experiments
% Input:
%   obj (CPepNormalization)
%       Normalization processor instance
%   checked_pep_file_path (1 x 1 char/string)
%       Path of the file containing checked peptides

fprintf('Starting requantifying normalization peptides...\n');

% Process each experiment
for i_exp = 1:length(obj.experimentNames)
    experimentName = obj.experimentNames{i_exp};
    fprintf('Processing experiment %d/%d: %s\n', i_exp, ...
        length(obj.experimentNames), experimentName);

    requant_processor = CMSMSPepDeconv( ...
        'modify.ini', ...
        '', '', ...     % modifications
        fullfile(obj.spectra_dir, experimentName), ...
        obj.ms1_tolerance, ...
        0.02, ...
        obj.alpha, ...
        '', '', ...     % fasta_path, regexp
        '', ...         % pepSpec_path
        1, 2, 0.5, 0.1, ...
        {}, ...         % enzyme
        fullfile(obj.result_dir, experimentName), ...
        [1,2], ...      % ionTypes
        checked_pep_file_path ...
        );
    requant_processor.requant_norm_pep();
end

fprintf('Requantification of normalization peptides completed.\n');
end