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
