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

% Input validation
if length(cstrIMP) ~= length(lfMasses) || length(cstrIMP) ~= length(abundance)
    error('CIMPGatherQuant:InvalidInput', ...
        'cstrIMP, lfMasses, and abundance must have the same length.');
end

% record the raw file if haven't been record
[obj, idx_raw] = obj.ensure_raw(raw_name);
raw = obj.get_raw(idx_raw);
raw.length = raw.length + 1;

% reallocate memory if there is no available place to record this
if raw.length > raw.capacity
    raw.curRts(raw.capacity+obj.m_buff_length, 1) = 0;
    raw.curIntens(raw.capacity+obj.m_buff_length, 1) = 0;
    raw.curMz(raw.capacity+obj.m_buff_length, 1) = 0;
    raw.curCharge(raw.capacity+obj.m_buff_length, 1) = 0;
    raw.ratioMatrix = [raw.ratioMatrix;zeros(obj.m_buff_length,size(raw.ratioMatrix,2))];
    raw.capacity = raw.capacity + obj.m_buff_length;
end

% record formally
raw.curRts(raw.length) = curRts;
raw.curIntens(raw.length) = curIntens;
raw.curMz(raw.length) = curMz;
raw.curCharge(raw.length) = cur_ch;
for iIso = 1:length(cstrIMP)
    if isKey(raw.mapIMPNames,cstrIMP{iIso})
        raw.ratioMatrix(raw.length,raw.mapIMPNames(cstrIMP{iIso})) = abundance(iIso);
    else
        % If not found, record it and append a column to ratioMatrix
        raw.mapIMPNames(cstrIMP{iIso}) = raw.mapIMPNames.Count+1;
        raw.impMass = [raw.impMass, lfMasses(iIso)];
        raw.impNames{raw.mapIMPNames.Count,1} = cstrIMP{iIso};
        raw.ratioMatrix = [raw.ratioMatrix,zeros(raw.capacity,1)];
        raw.ratioMatrix(raw.length,raw.mapIMPNames.Count) = abundance(iIso);
    end
end
obj = obj.set_raw(idx_raw, raw);
end

