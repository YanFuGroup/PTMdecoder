classdef CMSMSLevelService < handle
    % Stage service for MSMS-level quantification.

    properties (Access = private)
        m_msms_cfg
    end

    methods
        function obj = CMSMSLevelService(msms_cfg)
            % Input:
            %   msms_cfg (CMSMSPepDeconvConfig)
            %       config for MSMS processor
            if nargin < 1 || isempty(msms_cfg)
                error('CMSMSLevelService:MissingConfig', ...
                    'msms_cfg must be provided.');
            end
            obj.m_msms_cfg = msms_cfg;
        end

        function run(obj)
            % Run MSMS-level stage.
            processor = CMSMSPepDeconv(obj.m_msms_cfg);
            processor.runMsmsLevel();
        end
    end
end
