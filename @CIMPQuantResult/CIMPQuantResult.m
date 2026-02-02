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

        function rec = toStruct(obj)
            rec = struct(...
                'imp_name', obj.imp_name, ...
                'charge', obj.charge, ...
                'raw_name', obj.raw_name, ...
                'mass_center', obj.mass_center, ...
                'low_mz_bound', obj.low_mz_bound, ...
                'high_mz_bound', obj.high_mz_bound, ...
                'area', obj.area, ...
                'rt_peaks', obj.rt_peaks);
        end
    end

    methods (Static)
        function recs = toStructArray(results)
            if isempty(results)
                recs = repmat(struct('imp_name','', 'charge',0, 'raw_name','', ...
                    'mass_center',0, 'low_mz_bound',0, 'high_mz_bound',0, 'area',0, ...
                    'rt_peaks',struct('start',{},'end',{},'ratio',{},'check_label',{})), 0, 1);
                return;
            end
            recs = arrayfun(@(x) x.toStruct(), results);
        end

        function obj = empty(varargin)
            obj = CIMPQuantResult;
            obj = obj([]);
        end
    end
end
