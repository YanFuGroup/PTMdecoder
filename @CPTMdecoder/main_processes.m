function main_processes(~, varargin)
% Load the files and run PTMDecoder iteratively
% Input:
%   varargin (1 x N cell)
%       list of task parameter file paths

for i = 1:length(varargin)
    workflow_config = CPTMdecoderWorkflowConfig.fromParamFile(varargin{i});
    CLogger.info('Processing parameter file: %s', varargin{i});
    workflow_runner = CPTMdecoderWorkflowRunner(workflow_config);
    workflow_runner.run();
end
end