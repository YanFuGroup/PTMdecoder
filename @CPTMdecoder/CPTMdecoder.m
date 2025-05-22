classdef CPTMdecoder
    % Main workflow of CPTMDecoder

    methods
        function obj = CPTMdecoder(varargin)
            if nargin == 0
                % Default constructor
                error('No input arguments');

                % TODO: Add the default parameter file
                
            else
                % The arguments are all file paths
                % Check if the files exist
                obj.check_files_existence(varargin{:});

                obj.main_processes(varargin{:});
            end
        end

        % Check if the files exist
        check_files_existence(obj, varargin);

        % Load the files and run PTMDecoder iteratively
        main_processes(obj, varargin);
    end
end