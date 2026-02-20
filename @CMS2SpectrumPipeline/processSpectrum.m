function [bSuccess,cstrIMP,abundance,ionTypePosCharge,ionIntens,frageff,warning_msg,is_X_not_full_column_rank] = processSpectrum(obj)
% processSpectrum - Layered orchestration with rollback to legacy path
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

ionTypePosCharge = [];
ionIntens = [];
frageff = [];
is_X_not_full_column_rank = false;

try
    % Read spectrum and precursor information
    [obj.m_expPeaks,obj.m_iCharge,dPrecursorMZ] = obj.m_cMgfDatasetIO.read_oneSpec( ...
        obj.m_strDatasetName,obj.m_strSpecName);
    obj.m_dPrecursorMass = (dPrecursorMZ - CConstant.pmass) * obj.m_iCharge;

    % Build calculation context for mass-side modules
    mass_ctx = struct( ...
        'm_pepSeq', obj.m_pepSeq, ...
        'm_isProtN', obj.m_isProtN, ...
        'm_isProtC', obj.m_isProtC, ...
        'm_fixedModNameMass', {obj.m_fixedModNameMass}, ...
        'm_variableModNameMass', {obj.m_variableModNameMass}, ...
        'm_ms1_tolerance', obj.m_ms1_tolerance, ...
        'm_enzyme', obj.m_enzyme, ...
        'm_dPrecursorMass', obj.m_dPrecursorMass, ...
        'm_strSpecName', obj.m_strSpecName, ...
        'm_iCharge', obj.m_iCharge, ...
        'm_ionTypes', obj.m_ionTypes);

    % Enumerate feasible peptidoforms and corresponding non-redundant ions
    fixedPosMod = CMS2MassCalculator.getFixedPosMod(mass_ctx);
    [bSuccess,inxSites,massArrangement,warning_msg] = CMS2MassCalculator.getMassArrangement(mass_ctx,fixedPosMod);
    if ~bSuccess
        cstrIMP = [];
        abundance = [];
        return;
    end

    vNonRedunTheoryIonMz = CMS2MassCalculator.getNonRedunIons(mass_ctx,inxSites,massArrangement,fixedPosMod);

    % Match and preprocess peaks on non-redundant ion space
    matchedExpPeaks = CMS2PeakMatcher.match(obj.m_expPeaks,vNonRedunTheoryIonMz,obj.m_ms2_tolerance);
    matchedExpPeaks = CMS2PeakMatcher.processPeaks(matchedExpPeaks,obj.m_alpha);

    if isempty(matchedExpPeaks)
        if size(massArrangement,1) ~= 1
            cstrIMP = [];
            abundance = [];
            bSuccess = false;
            warning_msg = ['There is no non-redundant peak for discriminating the' ...
                ' peptidoforms for ',obj.m_pepSeq, ' in ', obj.m_strSpecName, '!\n'];
            return;
        else
            abundance = 1;
            cstrIMP = CMS2ResultIO.formatImpStrings(massArrangement,fixedPosMod,obj.m_variableModNameMass,inxSites,obj.m_pepSeq);
            return;
        end
    end

    % Post-match peptidoform filtering (remove empty / strict-subset candidates)
    [massArrangement, vNonRedunTheoryIonMz] = CMS2PeakMatcher.postMatchFilterPeptidoforms( ...
        matchedExpPeaks, massArrangement, vNonRedunTheoryIonMz);
    if size(massArrangement,1) == 1
        abundance = 1;
        cstrIMP = CMS2ResultIO.formatImpStrings(massArrangement,fixedPosMod,obj.m_variableModNameMass,inxSites,obj.m_pepSeq);
        bSuccess = true;
        return;
    end

    % Solve abundance according to configured model/method
    solver_cfg = struct( ...
        'model', obj.m_model, ...
        'method', obj.m_method, ...
        'lambda', obj.m_lambda, ...
        'case_penalty_intens', obj.m_case_penalty_intens, ...
        'grid_penalty_intens', obj.m_grid_penalty_intens, ...
        'case_OLS_intens_weight', obj.m_case_OLS_intens_weight);

    [abundance, frageff, ionTypePosCharge, ionIntens, is_X_not_full_column_rank] = ...
        CMS2QuantSolver.solve(vNonRedunTheoryIonMz, matchedExpPeaks, massArrangement, solver_cfg);

    % Final thresholding + normalization
    abundance(abundance<obj.m_resFilterThres*max(abundance))=0;
    abundance=abundance/(sum(abundance)+eps);

    cstrIMP = CMS2ResultIO.formatImpStrings(massArrangement,fixedPosMod,obj.m_variableModNameMass,inxSites,obj.m_pepSeq);
    bSuccess = true;
catch ME
    % Rollback path: fallback to legacy one-shot implementation
    legacy = CEachSpectrumLocQuant(obj.m_pepSeq,obj.m_isProtN,obj.m_isProtC, ...
        obj.m_cMgfDatasetIO,obj.m_strDatasetName,obj.m_strSpecName,obj.m_fixedModNameMass, ...
        obj.m_variableModNameMass,obj.m_model,obj.m_method,obj.m_lambda, ...
        obj.m_ms1_tolerance,obj.m_ms2_tolerance,obj.m_alpha,obj.m_resFilterThres,...
        obj.m_ionTypes,obj.m_enzyme,obj.m_case_penalty_intens,obj.m_grid_penalty_intens, ...
        obj.m_case_OLS_intens_weight);
    [bSuccess,cstrIMP,abundance,ionTypePosCharge,ionIntens,frageff,warning_msg,is_X_not_full_column_rank] = legacy.runEach();
    rollback_msg = ['[CMS2 rollback] ', ME.identifier, ': ', ME.message, '\n'];
    warning_msg = [warning_msg, rollback_msg];
end
end
