classdef CIMPRawIdentManager < handle
    % Manage per-raw identification stores for one peptide

    properties (Access=private)
        m_buff_length
        m_mapRawNames
        m_rawIdentData
    end

    methods
        function obj = CIMPRawIdentManager(buff_length)
            if nargin < 1
                buff_length = 50;
            end
            obj.m_buff_length = buff_length;
            obj.m_mapRawNames = containers.Map('KeyType','char','ValueType','any');
            obj.m_rawIdentData = {};
        end

        function rawStore = getOrCreate(obj, raw_name)
            idx_raw = obj.ensure_raw(raw_name);
            rawStore = obj.get_raw(idx_raw);
        end

        function setStore(obj, raw_name, rawStore)
            idx_raw = obj.ensure_raw(raw_name);
            obj.set_raw(idx_raw, rawStore);
        end

        function [raw_names, raw_ident_stores] = getEntries(obj)
            raw_names = obj.m_mapRawNames.keys;
            raw_ident_stores = cell(size(raw_names));
            for idx_keys = 1:obj.m_mapRawNames.Count
                idx_r = obj.m_mapRawNames(raw_names{idx_keys});
                raw_ident_stores{idx_keys} = obj.get_raw(idx_r);
            end
        end

        function clear(obj)
            obj.m_mapRawNames = containers.Map('KeyType','char','ValueType','any');
            obj.m_rawIdentData = {};
        end
    end

    methods (Access=private)
        function idx_raw = ensure_raw(obj, raw_name)
            if ~obj.m_mapRawNames.isKey(raw_name)
                idx_raw = obj.m_mapRawNames.Count + 1;
                obj.m_mapRawNames(raw_name) = idx_raw;
                obj.m_rawIdentData{idx_raw} = obj.init_raw();
            else
                idx_raw = obj.m_mapRawNames(raw_name);
            end
        end

        function raw = get_raw(obj, idx_raw)
            raw = obj.m_rawIdentData{idx_raw};
        end

        function set_raw(obj, idx_raw, raw)
            obj.m_rawIdentData{idx_raw} = raw;
        end

        function raw = init_raw(obj)
            raw = CIMPRawIdentStore(obj.m_buff_length);
        end
    end
end
