function obj = appendOneSpecQuant(obj,raw_name,curRts,curIntens,curMz,cur_ch,cstrIMP,lfMasses,abundance)
% Append one quantification result
% Input:
%   raw_name (1 x 1 char/string)
%       name of the run file (mgf)
%   curRts (1 x 1 double) minutes
%       retention time of this spectrum
%   curIntens (1 x 1 double) intensity
%       intensity of this spectrum
%   curMz (1 x 1 double) m/z
%       mass-to-charge ratio of this spectrum
%   cur_ch (1 x 1 double/int)
%       charge of the spectrum
%   cstrIMP (K x 1 cellstr/string)
%       string form of modified peptides (each is a different IMP)
%   lfMasses (1 x K double) or (K x 1 double) Da
%       masses of the modified peptides
%   abundance (K x 1 double)
%       relative abundance of each IMP for this spectrum

% record the raw file if haven't been record
if ~obj.m_mapRawNames.isKey(raw_name)
    idx_raw = obj.m_mapRawNames.Count + 1;
    obj.m_mapRawNames(raw_name) = idx_raw;
    obj.m_length{idx_raw} = 0;
    obj.m_capacity{idx_raw} = obj.m_buff_length;
    obj.m_curRts{idx_raw} = zeros(obj.m_buff_length,1);
    obj.m_curIntens{idx_raw} = zeros(obj.m_buff_length,1);
    obj.m_curMz{idx_raw} = zeros(obj.m_buff_length,1);
    obj.m_curCharge{idx_raw} = zeros(obj.m_buff_length,1);
    obj.m_ratioMatrix{idx_raw} = zeros(obj.m_buff_length,0);
    obj.m_cstrIMPNames{idx_raw} = cell(0);
    obj.m_IMPMass{idx_raw} = [];
    obj.m_mapIMPNames{idx_raw} = containers.Map();
end

% find the index of the raw file
idx_raw = obj.m_mapRawNames(raw_name);
obj.m_length{idx_raw} = obj.m_length{idx_raw}+1;

% reallocate memory if there is no available place to record this
if obj.m_length{idx_raw} > obj.m_capacity{idx_raw}
    obj.m_curRts{idx_raw}(obj.m_capacity{idx_raw}+obj.m_buff_length, 1) = 0;
    obj.m_curIntens{idx_raw}(obj.m_capacity{idx_raw}+obj.m_buff_length, 1) = 0;
    obj.m_curMz{idx_raw}(obj.m_capacity{idx_raw}+obj.m_buff_length, 1) = 0;
    obj.m_curCharge{idx_raw}(obj.m_capacity{idx_raw}+obj.m_buff_length, 1) = 0;
    obj.m_ratioMatrix{idx_raw} = [obj.m_ratioMatrix{idx_raw};zeros(obj.m_buff_length,size(obj.m_ratioMatrix{idx_raw},2))];
    obj.m_capacity{idx_raw} = obj.m_capacity{idx_raw}+obj.m_buff_length;
end

% record formally
obj.m_curRts{idx_raw}(obj.m_length{idx_raw}) = curRts;
obj.m_curIntens{idx_raw}(obj.m_length{idx_raw}) = curIntens;
obj.m_curMz{idx_raw}(obj.m_length{idx_raw}) = curMz;
obj.m_curCharge{idx_raw}(obj.m_length{idx_raw}) = cur_ch;
for iIso = 1:length(cstrIMP)
    if isKey(obj.m_mapIMPNames{idx_raw},cstrIMP{iIso})
        obj.m_ratioMatrix{idx_raw}(obj.m_length{idx_raw},...
            obj.m_mapIMPNames{idx_raw}(cstrIMP{iIso})) = abundance(iIso);
    else
        % If not found, record it and append a column to ratioMatrix
        obj.m_mapIMPNames{idx_raw}(cstrIMP{iIso}) = obj.m_mapIMPNames{idx_raw}.Count+1;
        obj.m_IMPMass{idx_raw} = [obj.m_IMPMass{idx_raw}, lfMasses(iIso)];
        obj.m_cstrIMPNames{idx_raw}{obj.m_mapIMPNames{idx_raw}.Count,1} = cstrIMP{iIso};
        obj.m_ratioMatrix{idx_raw} = [obj.m_ratioMatrix{idx_raw},zeros(obj.m_capacity{idx_raw},1)];
        obj.m_ratioMatrix{idx_raw}(obj.m_length{idx_raw},obj.m_mapIMPNames{idx_raw}.Count) = abundance(iIso);
    end
end
end

