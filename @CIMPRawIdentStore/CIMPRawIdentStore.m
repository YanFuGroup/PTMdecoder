classdef CIMPRawIdentStore
    % Raw identification store for one raw file
    % Holds per-spectrum data and IMP ratio matrix with dynamic growth

    properties
        length
        capacity
        curRts
        curIntens
        curMz
        curCharge
        ratioMatrix
        impNames
        impMass
        mapIMPNames
        buffLength
    end

    methods
        function obj = CIMPRawIdentStore(buffLength)
            if nargin < 1
                buffLength = 50;
            end
            obj.buffLength = buffLength;
            obj.length = 0;
            obj.capacity = buffLength;
            obj.curRts = zeros(buffLength,1);
            obj.curIntens = zeros(buffLength,1);
            obj.curMz = zeros(buffLength,1);
            obj.curCharge = zeros(buffLength,1);
            obj.ratioMatrix = zeros(buffLength,0);
            obj.impNames = cell(0,1);
            obj.impMass = [];
            obj.mapIMPNames = containers.Map();
        end

        function obj = appendSpecQuant(obj, curRts, curIntens, curMz, cur_ch, cstrIMP, lfMasses, abundance)
            if numel(cstrIMP) ~= numel(lfMasses) || numel(cstrIMP) ~= numel(abundance)
                error('CIMPRawIdentStore:InvalidInput', ...
                    'cstrIMP, lfMasses, and abundance must have the same length.');
            end

            obj.length = obj.length + 1;
            if obj.length > obj.capacity
                obj.curRts(obj.capacity + obj.buffLength, 1) = 0;
                obj.curIntens(obj.capacity + obj.buffLength, 1) = 0;
                obj.curMz(obj.capacity + obj.buffLength, 1) = 0;
                obj.curCharge(obj.capacity + obj.buffLength, 1) = 0;
                obj.ratioMatrix = [obj.ratioMatrix; zeros(obj.buffLength, size(obj.ratioMatrix,2))];
                obj.capacity = obj.capacity + obj.buffLength;
            end

            obj.curRts(obj.length) = curRts;
            obj.curIntens(obj.length) = curIntens;
            obj.curMz(obj.length) = curMz;
            obj.curCharge(obj.length) = cur_ch;

            for idx = 1:numel(cstrIMP)
                imp_name = cstrIMP{idx};
                if isKey(obj.mapIMPNames, imp_name)
                    col_idx = obj.mapIMPNames(imp_name);
                    obj.ratioMatrix(obj.length, col_idx) = abundance(idx);
                else
                    new_idx = obj.mapIMPNames.Count + 1;
                    obj.mapIMPNames(imp_name) = new_idx;
                    obj.impMass = [obj.impMass, lfMasses(idx)];
                    obj.impNames{new_idx,1} = imp_name;
                    obj.ratioMatrix = [obj.ratioMatrix, zeros(obj.capacity,1)];
                    obj.ratioMatrix(obj.length, new_idx) = abundance(idx);
                end
            end
        end
    end
end
