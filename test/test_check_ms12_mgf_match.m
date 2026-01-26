classdef test_check_ms12_mgf_match < matlab.unittest.TestCase
    % Test checking MS1/MGF filename matching logic in CMSMSPepDeconv
    
    properties
        RepoRoot
        TestDir
        Obj
    end
    
    methods(TestClassSetup)
        function setupEnvironment(testCase)
            % Setup paths
            testCase.RepoRoot = fileparts(fileparts(mfilename('fullpath'))); % project dir
            addpath(genpath(testCase.RepoRoot)); % Add all code to path
            
            % Create temporary directory for test files
            testCase.TestDir = fullfile(testCase.RepoRoot, 'test', 'temp_ms12_test');
            if exist(testCase.TestDir, 'dir')
                rmdir(testCase.TestDir, 's');
            end
            mkdir(testCase.TestDir);
            
            % Create Object
            mockParam = testCase.createMockTaskParam();
            % Instantiate the class
            testCase.Obj = CMSMSPepDeconv(mockParam);
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
            testCase.verifyWarningFree(@() testCase.Obj.check_whether_ms12_mgf_name_match());
        end
        
        function testSuffixMatchHCDFT(testCase)
            % Test 2: Suffix Match (HCDFT)
            % file2_HCDFT.mgf <-> file2.ms1
            testCase.createFile('file2_HCDFT.mgf');
            testCase.createFile('file2.ms1');
            
            % Function should run without error
            testCase.verifyWarningFree(@() testCase.Obj.check_whether_ms12_mgf_name_match());
        end
        
        function testSuffixMatchCIDIT(testCase)
            % Test 3: Suffix Match (CIDIT)
            % file3_CIDIT.mgf <-> file3.ms1
            testCase.createFile('file3_CIDIT.mgf');
            testCase.createFile('file3.ms1');
            
            % Function should run without error
            testCase.verifyWarningFree(@() testCase.Obj.check_whether_ms12_mgf_name_match());
        end
        
        function testMissingMS1(testCase)
            % Test 4: Missing MS1 (Should Error)
            % file4.mgf -> no ms1
            testCase.createFile('file4.mgf');
            
            try
                testCase.Obj.check_whether_ms12_mgf_name_match();
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
            testCase.verifyWarningFree(@() testCase.Obj.check_whether_ms12_mgf_name_match());
        end

        function testReviewPSM_NoMatch(testCase)
            % Test 6: Suffix present but no matching MS1 file. (Should Error)
            testCase.createFile('file6_HCDFT.mgf');
            testCase.createFile('file5_HCDFT.ms1');
            
            try
                testCase.Obj.check_whether_ms12_mgf_name_match();
                testCase.verifyFail('Expected error was not thrown');
            catch ME
                testCase.verifyTrue(contains(ME.message, 'does not have a corresponding .ms1 file'), ...
                    'Error message did not contain expected text.');
            end
        end
    end
    
    methods (Access = private)
        function mockParam = createMockTaskParam(testCase)
            mockParam.m_spec_dir_path = testCase.TestDir;
            mockParam.m_mod_file_path = fullfile(testCase.RepoRoot, 'modify.ini'); % Point to existing file
            mockParam.m_fasta_file_path = '';
            mockParam.m_regular_express = '';
            mockParam.m_pep_spec_file_path = '';
            mockParam.m_filtered_res_file_path = '';
            mockParam.m_output_dir_path = '';
            mockParam.m_model = '';
            mockParam.m_method = '';
            mockParam.m_lambda = 0;
            mockParam.m_ms1_tolerance = 0;
            mockParam.m_ms2_tolerance = 0;
            mockParam.m_alpha = 0;
            mockParam.m_result_filter_threshold = 0;
            mockParam.m_enzyme_name = 'Trypsin';
            mockParam.m_enzyme_limits = [];
            mockParam.m_checked_peptides_res_path = '';
            mockParam.m_msms_res_path = '';
            mockParam.m_min_MSMS_num = 1;
            mockParam.m_fixed_mod = {};
            mockParam.m_variable_mod = {};
        end
        
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
