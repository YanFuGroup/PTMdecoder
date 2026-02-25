classdef CPeptideLevelPipeline < handle
    % Orchestration layer for peptide-level quantification workflows.
    % It delegates concrete execution to executor-layer components.

    properties (Access=private)
        m_processing_executor
        m_align_requant_executor
    end

    methods
        function obj = CPeptideLevelPipeline(processing_executor, align_requant_executor)
            if nargin < 1
                processing_executor = [];
            end
            if nargin < 2
                align_requant_executor = [];
            end
            obj.m_processing_executor = processing_executor;
            obj.m_align_requant_executor = align_requant_executor;
        end

        function setProcessingExecutor(obj, processing_executor)
            obj.m_processing_executor = processing_executor;
        end

        function setAlignRequantExecutor(obj, align_requant_executor)
            obj.m_align_requant_executor = align_requant_executor;
        end

        function block = quantifyPeptideBlock(obj, prot_names_pos, rawIdentManager)
            obj.assertProcessingExecutor();
            block = obj.m_processing_executor.quantifyPeptideBlock(prot_names_pos, rawIdentManager);
        end

        function block = requantifyPeptideBlock(obj, prot_names_pos, rawIdentManager, pep_rtrange_map)
            obj.assertProcessingExecutor();
            block = obj.m_processing_executor.requantifyPeptideBlock(prot_names_pos, rawIdentManager, pep_rtrange_map);
        end

        function drawImpXicForBlock(obj, rawIdentManager, pep_rtrange_map, dir_save, color_map, legend_map)
            obj.assertProcessingExecutor();
            obj.m_processing_executor.drawImpXicForBlock(rawIdentManager, pep_rtrange_map, dir_save, color_map, legend_map);
        end

        function [pep_rtrange_map, report] = buildAlignedRtRangeMap(obj, msms_result, fdr_filtered_result_path, buildRawIdentManagerFn)
            obj.assertAlignRequantExecutor();
            [pep_rtrange_map, report] = obj.m_align_requant_executor.buildAlignedRtRangeMap(msms_result, fdr_filtered_result_path, buildRawIdentManagerFn);
        end

        function writeAlignmentReport(obj, report, output_path)
            obj.assertAlignRequantExecutor();
            obj.m_align_requant_executor.writeAlignmentReport(report, output_path);
        end
    end

    methods (Access=private)
        function assertProcessingExecutor(obj)
            if isempty(obj.m_processing_executor)
                error('CPeptideLevelPipeline:MissingProcessingExecutor', ...
                    'processing executor must be provided before quant/requant/draw.');
            end
        end

        function assertAlignRequantExecutor(obj)
            if isempty(obj.m_align_requant_executor)
                error('CPeptideLevelPipeline:MissingAlignRequantExecutor', ...
                    'align-requant executor must be provided before alignment steps.');
            end
        end
    end
end
