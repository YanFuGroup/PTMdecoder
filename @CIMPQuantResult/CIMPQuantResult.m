classdef CIMPQuantResult
    % Data model for quantification result of an IMP

    properties
        imp_name
        charge
        raw_name
        mass_center
        low_mz_bound
        high_mz_bound
        area
        rt_peaks
    end

    methods
        function obj = CIMPQuantResult(imp_name, charge, raw_name, mass_center, low_mz_bound, high_mz_bound, area, rt_peaks)
            % Constructor for CIMPQuantResult
            % input:
            %   imp_name (1 x 1 char/string)
            %       the name of the IMP
            %   charge (1 x 1 double/int)
            %       charge state of the IMP
            %   raw_name (1 x 1 char/string)
            %       the name of the raw data file
            %   mass_center (1 x 1 double) m/z
            %       central mass value
            %   low_mz_bound (1 x 1 double) m/z
            %       low precursor m/z bound
            %   high_mz_bound (1 x 1 double) m/z
            %       high precursor m/z bound
            %   area (1 x 1 double) area
            %       quantified area
            %   rt_peaks (struct array) 
            %       array of structs, each with fields:
            %       - rt_start (double) minutes: start retention time of the peak
            %       - rt_end (double) minutes: end retention time of the peak
            %       - ratio (double): ratio contribution of the peak
            %       - check_label (int): 1 if selected peak, 0 otherwise,and 2 for checked selected peak
            if nargin == 0
                return;
            end
            obj.imp_name = imp_name;
            obj.charge = charge;
            obj.raw_name = raw_name;
            obj.mass_center = mass_center;
            obj.low_mz_bound = low_mz_bound;
            obj.high_mz_bound = high_mz_bound;
            obj.area = area;
            obj.rt_peaks = rt_peaks;
        end
    end

    methods (Static)
        function obj = empty(varargin)
            obj = CIMPQuantResult;
            obj = obj([]);
        end
    end
end
