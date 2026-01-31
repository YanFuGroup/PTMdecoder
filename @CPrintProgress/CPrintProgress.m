classdef CPrintProgress
    % For printing the progress of the long processing
    
    properties
        m_total_step;
        m_progress_1_100 = 1; % record the current progress
        m_string_length = 0;% record the last printing length
    end
    
    methods
        function obj = CPrintProgress(total_step)
            % Input:
            %   total_step (1 x 1 double/int)
            %       total number of steps
            obj.m_total_step = total_step;
        end
        
        function obj = update_show(obj,current_step)
            % Update the current status and show the progress
            % Input:
            %   current_step (1 x 1 double/int)
            %       current step number
            if double(current_step)/double(obj.m_total_step) > obj.m_progress_1_100/100
                progress_string = repmat('\b',1,obj.m_string_length);
                fprintf(progress_string);
                obj.m_string_length = fprintf('%d%%',obj.m_progress_1_100);
                obj.m_progress_1_100 = ceil(double(current_step)/double(obj.m_total_step)*100);
            end
        end

        function last_update(obj)
            % Update the last time
            progress_string = repmat('\b',1,obj.m_string_length);
            fprintf(progress_string);
        end
    end
end

