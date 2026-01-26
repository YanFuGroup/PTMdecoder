function system_test_runner()
    fclose all;

    % Set path
    currentDir = fileparts(mfilename('fullpath'));
    projectDir = fileparts(currentDir); % Go back one level to the code root directory
    addpath(projectDir);
    
    % Configure test data
    testDataDir = fullfile(currentDir, 'data');
    msms_pep_paramFile = fullfile(testDataDir, 'msms_pep_site.param');
    pep_requant_paramFile = fullfile(testDataDir, 'requant_pep_site.param');
    goldenDir = fullfile(currentDir, 'golden');
    outputDir = fullfile(currentDir, 'output');
    
    % Clean up old output (optional, to prevent interference)
    if isfolder(outputDir)
        try
            rmdir(outputDir, 's');
        catch
            fprintf('Warning: Could not remove output directory. Attempting to delete files inside.\n');
            delete(fullfile(outputDir, '*'));
        end
    end
    
    if ~isfolder(outputDir)
        mkdir(outputDir);
    end
    
    fprintf('=== Start running tests ===\n');
    
    fprintf('Running msms-pep procedure...\n');
    main(msms_pep_paramFile); 

    fprintf('Running pep-requant procedure...\n');
    main(pep_requant_paramFile);
    
    fprintf('Running drawXIC test...\n');
    helper_test_draw_xic(projectDir, testDataDir, outputDir);

    fprintf('Run completed, start comparing results...\n');
    
    % Comparison section
    allMatch = compare_dir_recursive(goldenDir, goldenDir, outputDir);
    
    if allMatch
        fprintf('\n=== Test passed: Refactoring safe ===\n');
    else
        fprintf('\n=== Test failed: Please check code changes ===\n');
    end
end

function isMatch = compare_binary_files(file1, file2)
    % Read and compare the contents of two files byte by byte
    f1 = java.io.File(file1);
    f2 = java.io.File(file2);
    
    if f1.length() ~= f2.length()
        isMatch = false; return;
    end
    
    fid1 = fopen(file1, 'r');
    content1 = fread(fid1, inf, 'uint8');
    fclose(fid1);
    
    fid2 = fopen(file2, 'r');
    content2 = fread(fid2, inf, 'uint8');
    fclose(fid2);
    
    isMatch = isequal(content1, content2);
end

function isMatch = compare_text_files(file1, file2, tol)
    % Compare text files allowing for small numerical differences
    
    % Read files
    str1 = fileread(file1);
    str2 = fileread(file2);
    
    % If identical strings, pass immediately
    if strcmp(str1, str2)
        isMatch = true;
        return;
    end
    
    % 1. Extract structural parts (non-numbers) to ensure template matches
    % Regex to separate numbers from text
    % Pattern matches ints, floats, sci-notation: -1.23e-4, .5, 100
    numPattern = '[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?';
    
    [nums1, split1] = regexp(str1, numPattern, 'match', 'split');
    [nums2, split2] = regexp(str2, numPattern, 'match', 'split');
    
    % Compare non-numeric structure (the headers, labels, tags)
    if length(split1) ~= length(split2) || ~isequal(split1, split2)
        % Structure differs
        isMatch = false;
        return;
    end
    
    % 2. Compare the extracted numbers with tolerance
    n1 = str2double(nums1);
    n2 = str2double(nums2);
    
    % Filter out NaNs if any (str2double might produce NaN for weird tokens if regex is loose, 
    % but current regex is strict enough for numbers. 
    % However, strict inequality check handles NaNs (NaN ~= NaN). 
    % We should treat NaNs as equal if both are NaN.)
    
    if length(n1) ~= length(n2)
        isMatch = false;
        return;
    end
    
    diffs = abs(n1 - n2);
    
    % Check tolerance
    % Handle NaNs: if both are NaN, diff is usually NaN. 
    % We need to ensure we compare values where not both are NaN.
    
    validMask = ~isnan(n1) & ~isnan(n2);
    
    % If one is NaN and other is not, it's a mismatch
    if any(isnan(n1) ~= isnan(n2))
        isMatch = false;
        return;
    end
    
    if any(diffs(validMask) > tol)
        isMatch = false;
        return;
    end
    
    isMatch = true;
end

function allMatch = compare_dir_recursive(currentDir, rootGoldenDir, rootOutputDir)
    files = dir(currentDir);
    allMatch = true;
    
    for i = 1:length(files)
        if strcmp(files(i).name, '.') || strcmp(files(i).name, '..')
            continue;
        end
        
        fullGoldenPath = fullfile(currentDir, files(i).name);
        % Calculate relative path assuming rootGoldenDir does not end with separator
        % and fullGoldenPath starts with rootGoldenDir + separator
        relativePath = fullGoldenPath(length(rootGoldenDir)+2:end);
        
        if files(i).isdir
            if ~compare_dir_recursive(fullGoldenPath, rootGoldenDir, rootOutputDir)
                allMatch = false;
            end
        else
            outputPath = fullfile(rootOutputDir, relativePath);
            
            if ~isfile(outputPath)
                fprintf('[Missing] New result did not generate file: %s\n', relativePath);
                allMatch = false;
                continue;
            end
            
            [~, ~, ext] = fileparts(files(i).name);
            isTextFile = any(strcmpi(ext, {'.txt'}));
            
            if isTextFile
                passed = compare_text_files(fullGoldenPath, outputPath, 1e-5);
            else
                passed = compare_binary_files(fullGoldenPath, outputPath);
            end
    
            if passed
                fprintf('[Pass] %s\n', relativePath);
            else
                fprintf('[Fail] %s content inconsistent!\n', relativePath);
                allMatch = false;
            end
        end
    end
end