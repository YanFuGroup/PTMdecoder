classdef CPeptideQuantService < handle
    % Stage service for peptide-level quantification.

    properties (Access = private)
        m_msms_cfg
    end

    methods
        function obj = CPeptideQuantService(msms_cfg)
            % Input:
            %   msms_cfg (CMSMSPepDeconvConfig)
            %       config for peptide-level processor
            if nargin < 1 || isempty(msms_cfg)
                error('CPeptideQuantService:MissingConfig', ...
                    'msms_cfg must be provided.');
            end
            obj.m_msms_cfg = msms_cfg;
        end

        function run(obj)
            % Run peptide-level quantification stage.
            processor = CMSMSPepDeconv(obj.m_msms_cfg);
            processor.runImpQuantLevel();
        end
    end
end
