function obj=delUselessRaws(obj)
% Delete useless raws
% Input:
%   obj (CIMPGatherQuant)
%       Quantification aggregator instance
% Output:
%   obj (CIMPGatherQuant)
%       Updated instance with trimmed buffers
for idx = 1:length(obj.m_curRts)
    obj.m_curRts{idx}(obj.m_length{idx}+1:end,:)=[];
    obj.m_curIntens{idx}(obj.m_length{idx}+1:end,:)=[];
    obj.m_curMz{idx}(obj.m_length{idx}+1:end,:)=[];
    obj.m_curCharge{idx}(obj.m_length{idx}+1:end,:)=[];
    obj.m_ratioMatrix{idx}(obj.m_length{idx}+1:end,:)=[];
end
end

