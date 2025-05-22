function obj = appendOneSpecQuant(obj,raw_name,curRts,curIntens,curMz,cur_ch,cstrPepIso,lfMasses,abundance)
% Append one quantification result
% Input: raw_name is the name of the run file (mgf)
%       curRts is the retention time of this spectrum
%       curIntens is the intensity of this spectrum
%       curMz is the mass-to-charge ratio of this spectrum
%       cur_ch is the charge of the spectrum
%       cstrPepIso is the string form of the modified peptide, a column vector, each is a different modified peptide
%       lfMass is the mass of the modified peptide
%       abundance is the result of single spectrum quantification, a column vector, each represents the relative abundance of the above modified peptide

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
for iIso = 1:length(cstrPepIso)
    if isKey(obj.m_mapIMPNames{idx_raw},cstrPepIso{iIso})
        obj.m_ratioMatrix{idx_raw}(obj.m_length{idx_raw},...
            obj.m_mapIMPNames{idx_raw}(cstrPepIso{iIso})) = abundance(iIso);
    else
        % If not found, record it and append a column to ratioMatrix
        obj.m_mapIMPNames{idx_raw}(cstrPepIso{iIso}) = obj.m_mapIMPNames{idx_raw}.Count+1;
        obj.m_IMPMass{idx_raw} = [obj.m_IMPMass{idx_raw}, lfMasses(iIso)];
        obj.m_cstrIMPNames{idx_raw}{obj.m_mapIMPNames{idx_raw}.Count,1} = cstrPepIso{iIso};
        obj.m_ratioMatrix{idx_raw} = [obj.m_ratioMatrix{idx_raw},zeros(obj.m_capacity{idx_raw},1)];
        obj.m_ratioMatrix{idx_raw}(obj.m_length{idx_raw},obj.m_mapIMPNames{idx_raw}.Count) = abundance(iIso);
    end
end
end

