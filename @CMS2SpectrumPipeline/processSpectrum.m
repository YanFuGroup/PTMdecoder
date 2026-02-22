function [bSuccess,cstrIMP,abundance,ionTypePosCharge,ionIntens,frageff,warning_msg,is_X_not_full_column_rank] = processSpectrum(obj)
% processSpectrum - Layered orchestration for MS2 spectrum quantification
% Outputs:
%   bSuccess (1 x 1 logical)
%       Whether processing succeeded.
%   cstrIMP (K x 1 cell)
%       IMP identifiers in string form.
%   abundance (K x 1 double)
%       Relative abundance of IMPs.
%   ionTypePosCharge (L x 3 double)
%       Fragment ion info [type, position, charge].
%   ionIntens (L x 1 double)
%       Intensity for ions in ionTypePosCharge.
%   frageff (L x 1 double)
%       Fragmentation efficiency for ions in ionTypePosCharge.
%   warning_msg (1 x 1 char/string)
%       Warning text generated during processing.
%   is_X_not_full_column_rank (1 x 1 logical)
%       True if X is rank-deficient; false otherwise.

try
    % Read spectrum and precursor information
    [expPeaks,iCharge,precursorMZ] = obj.m_cMgfDatasetIO.read_oneSpec( ...
        obj.m_strDatasetName,obj.m_strSpecName);

    peptideCtx = struct( ...
        'pepSeq', obj.m_pepSeq, ...
        'isProtN', obj.m_isProtN, ...
        'isProtC', obj.m_isProtC, ...
        'fixedModNameMass', {obj.m_fixedModNameMass}, ...
        'variableModNameMass', {obj.m_variableModNameMass});

    spectrumCtx = struct( ...
        'datasetName', obj.m_strDatasetName, ...
        'specName', obj.m_strSpecName, ...
        'expPeaks', expPeaks, ...
        'iCharge', iCharge, ...
        'precursorMZ', precursorMZ);

    [bSuccess,cstrIMP,abundance,ionTypePosCharge,ionIntens,frageff,warning_msg,is_X_not_full_column_rank] = ...
        obj.processSpectrumWithContext(peptideCtx, spectrumCtx);
catch ME
    bSuccess = false;
    cstrIMP = [];
    abundance = [];
    ionTypePosCharge = [];
    ionIntens = [];
    frageff = [];
    is_X_not_full_column_rank = false;
    warning_msg = ['[CMS2] ', ME.identifier, ': ', ME.message, '\n'];
end
end
