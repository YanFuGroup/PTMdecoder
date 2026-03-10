classdef CPrintProgress
    % For printing the progress of the long processing
    
    properties
        m_total_step;
        m_progress_1_100 = 1; % record the next percentage threshold to log
        m_label = 'progress';
    end
    
    methods
        function obj = CPrintProgress(total_step, label)
            % Input:
            %   total_step (1 x 1 double/int)
            %       total number of steps
            %   label (1 x 1 char/string, optional)
            %       progress label written into logs
            obj.m_total_step = total_step;
            if nargin >= 2 && ~isempty(label)
                obj.m_label = char(label);
            end
        end
        
        function obj = update_show(obj,current_step)
            % Update the current status and show the progress
            % Input:
            %   current_step (1 x 1 double/int)
            %       current step number
            if obj.m_total_step <= 0
                return;
            end
            new_progress = floor(double(current_step) / double(obj.m_total_step) * 100);
            new_progress = min(max(new_progress, 0), 100);
            if new_progress >= obj.m_progress_1_100
                CLogger.progress(obj.m_label, current_step, obj.m_total_step);
                obj.m_progress_1_100 = new_progress + 1;
            end
        end

        function last_update(obj)
            % Update the last time
            if obj.m_total_step <= 0
                return;
            end
            % Avoid emitting duplicated 100%% when update_show already reached completion.
            if obj.m_progress_1_100 <= 100
                CLogger.progress(obj.m_label, obj.m_total_step, obj.m_total_step);
                obj.m_progress_1_100 = 101;
            end
        end
    end
end

