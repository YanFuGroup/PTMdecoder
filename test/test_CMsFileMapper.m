classdef test_CMsFileMapper < matlab.unittest.TestCase
    % Test checking MS1/MGF filename matching logic in CMsFileMapper
    
    properties
        RepoRoot
        TestDir
    end
    
    methods(TestClassSetup)
        function setupEnvironment(testCase)
            % Setup paths
            testCase.RepoRoot = fileparts(fileparts(mfilename('fullpath'))); % project dir
            addpath(genpath(testCase.RepoRoot)); % Add all code to path
            
            % Create temporary directory for test files
            testCase.TestDir = fullfile(testCase.RepoRoot, 'test', 'temp_ms12_test_mapper');
            if exist(testCase.TestDir, 'dir')
                rmdir(testCase.TestDir, 's');
            end
            mkdir(testCase.TestDir);
        end
    end
    
    methods(TestClassTeardown)
        function cleanupEnvironment(testCase)
            if exist(testCase.TestDir, 'dir')
                rmdir(testCase.TestDir, 's');
            end
        end
    end
    
    methods(TestMethodSetup)
        function resetTestDir(testCase)
            % Clear files in testDir before each test
            files = dir(fullfile(testCase.TestDir, '*.*'));
            for i = 1:length(files)
                if ~files(i).isdir
                    delete(fullfile(testCase.TestDir, files(i).name));
                end
            end
        end
    end
    
    methods(Test)
        function testDirectMatch(testCase)
            % Test 1: Direct Match
            % file1.mgf <-> file1.ms1
            testCase.createFile('file1.mgf');
            testCase.createFile('file1.ms1');
            
            % Function should run without error
            testCase.verifyWarningFree(@() CMsFileMapper(testCase.TestDir));
        end
        
        function testSuffixMatchHCDFT(testCase)
            % Test 2: Suffix Match (HCDFT)
            % file2_HCDFT.mgf <-> file2.ms1
            testCase.createFile('file2_HCDFT.mgf');
            testCase.createFile('file2.ms1');
            
            % Function should run without error
            testCase.verifyWarningFree(@() CMsFileMapper(testCase.TestDir));
        end
        
        function testSuffixMatchCIDIT(testCase)
            % Test 3: Suffix Match (CIDIT)
            % file3_CIDIT.mgf <-> file3.ms1
            testCase.createFile('file3_CIDIT.mgf');
            testCase.createFile('file3.ms1');
            
            % Function should run without error
            testCase.verifyWarningFree(@() CMsFileMapper(testCase.TestDir));
        end
        
        function testMissingMS1(testCase)
            % Test 4: Missing MS1 (Should Error)
            % file4.mgf -> no ms1
            testCase.createFile('file4.mgf');
            
            try
                CMsFileMapper(testCase.TestDir);
                testCase.verifyFail('Expected error was not thrown');
            catch ME
                testCase.verifyTrue(contains(ME.message, 'does not have a corresponding .ms1 file'), ...
                    'Error message did not contain expected text.');
            end
        end
        
        function testSuffixMismatchOnlyIfBaseMissing(testCase)
            % Test 5: Suffix present. Base MS1 missing, BUT full name MS1 exists. (Should PASS now)
            % file5_HCDFT.mgf -> file5_HCDFT.ms1 (Full match)
            % This was previously an error case, now it should be allowed.
            testCase.createFile('file5_HCDFT.mgf');
            testCase.createFile('file5_HCDFT.ms1'); 
            
            % Function should run without error
            testCase.verifyWarningFree(@() CMsFileMapper(testCase.TestDir));
        end

        function testReviewPSM_NoMatch(testCase)
            % Test 6: Suffix present but no matching MS1 file. (Should Error)
            testCase.createFile('file6_HCDFT.mgf');
            testCase.createFile('file5_HCDFT.ms1'); % Mismatch name
            
            try
                CMsFileMapper(testCase.TestDir);
                testCase.verifyFail('Expected error was not thrown');
            catch ME
                testCase.verifyTrue(contains(ME.message, 'does not have a corresponding .ms1 file'), ...
                    'Error message did not contain expected text.');
            end
        end
        
        function testMappingCorrectness(testCase)
           % Test 7: Verify that the mapping is actually correct
           testCase.createFile('file7_HCDFT.mgf');
           testCase.createFile('file7.ms1');
           
           mapper = CMsFileMapper(testCase.TestDir);
           ms1_name = mapper.get_ms1_stem('file7_HCDFT');
           testCase.verifyEqual(ms1_name, 'file7');
        end
    end
    
    methods (Access = private)
        function createFile(testCase, name)
            fullpath = fullfile(testCase.TestDir, name);
            fid = fopen(fullpath, 'w');
            if fid == -1
                error('Cannot create file %s', fullpath);
            end
            fprintf(fid, 'test content');
            fclose(fid);
        end
    end
end
