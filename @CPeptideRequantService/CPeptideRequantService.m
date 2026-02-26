classdef CPeptideRequantService < handle
    % Stage service for peptide-level re-quantification.

    properties (Access = private)
        m_msms_cfg
    end

    methods
        function obj = CPeptideRequantService(msms_cfg)
            % Input:
            %   msms_cfg (CMSMSPepDeconvConfig)
            %       config for peptide-level re-quant processor
            if nargin < 1 || isempty(msms_cfg)
                error('CPeptideRequantService:MissingConfig', ...
                    'msms_cfg must be provided.');
            end
            obj.m_msms_cfg = msms_cfg;
        end

        function run(obj)
            % Run peptide-level re-quantification stage.
            processor = CMSMSPepDeconv(obj.m_msms_cfg);
            processor.runImpRequantLevel();
        end
    end
end
