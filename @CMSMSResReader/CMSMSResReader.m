classdef CMSMSResReader
    % Reader for MSMS level result files
    % Usage:
    %   reader = CMSMSResReader();
    %   resultObj = reader.read_from_msms_res_file(msms_res_path);
    %   % resultObj is an instance of CMSMSResult

    
    methods
        function obj = CMSMSResReader()
        end

        % Read from a msms result file
        % Returns a CMSMSResult object
        resultObj = read_from_msms_res_file(obj, msms_res_path);

    end
end