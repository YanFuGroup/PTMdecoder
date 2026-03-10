function [bSuccess,cstrIMP,abundance,ionTypePosCharge,ionIntens,frageff,is_X_not_full_column_rank] = processSpectrumWithContext(obj, peptideCtx, spectrumCtx)
% processSpectrumWithContext - Process one spectrum by explicit input contexts
% Input:
%   obj (CMS2SpectrumPipeline)
%       Processor instance carrying static solver/config parameters.
%   peptideCtx (1 x 1 struct)
%       Peptide-side context for one processing call.
%       Required fields:
%           pepSeq (1 x N char/string)
%               Peptide sequence.
%           isProtN (1 x 1 logical)
%               Whether peptide is protein N-terminal.
%           isProtC (1 x 1 logical)
%               Whether peptide is protein C-terminal.
%           fixedModNameMass (1 x K cell)
%               Fixed modification definition.
%           variableModNameMass (1 x M cell)
%               Variable modification definition.
%   spectrumCtx (1 x 1 struct)
%       Spectrum-side context for one processing call.
%       Required fields:
%           datasetName (1 x N char/string)
%               Dataset(raw file) name of the spectrum.
%           specName (1 x N char/string)
%               Spectrum title/index name.
%           expPeaks (P x 2 double)
%               Experimental peaks in [mz, intensity] form.
%           iCharge (1 x 1 double)
%               Precursor charge state.
%           precursorMZ (1 x 1 double)
%               Precursor m/z.
% Output:
%   bSuccess (1 x 1 logical)
%       True if this spectrum is processed successfully; false otherwise.
%   cstrIMP (Q x 1 cell)
%       Candidate peptidoform strings.
%   abundance (Q x 1 double)
%       Relative abundance for each peptidoform in cstrIMP.
%   ionTypePosCharge (R x 3 double)
%       Ion descriptors [ionType, ionPos, ionCharge].
%   ionIntens (R x 1 double)
%       Ion intensity corresponding to ionTypePosCharge.
%   frageff (R x 1 double)
%       Fragment efficiency corresponding to ionTypePosCharge.
%   is_X_not_full_column_rank (1 x 1 logical)
%       Whether matrix X is rank-deficient during abundance solve.

ionTypePosCharge = [];
ionIntens = [];
frageff = [];
is_X_not_full_column_rank = false;

obj.m_pepSeq = peptideCtx.pepSeq;
obj.m_isProtN = peptideCtx.isProtN;
obj.m_isProtC = peptideCtx.isProtC;
obj.m_fixedModNameMass = peptideCtx.fixedModNameMass;
obj.m_variableModNameMass = peptideCtx.variableModNameMass;

obj.m_strDatasetName = spectrumCtx.datasetName;
obj.m_strSpecName = spectrumCtx.specName;
obj.m_expPeaks = spectrumCtx.expPeaks;
obj.m_iCharge = spectrumCtx.iCharge;
obj.m_dPrecursorMass = (spectrumCtx.precursorMZ - CConstant.pmass) * obj.m_iCharge;

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
[bSuccess,inxSites,massArrangement] = CMS2MassCalculator.getMassArrangement(mass_ctx,fixedPosMod);
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
        CLogger.debug(['There is no non-redundant peak for discriminating the' ...
            ' peptidoforms for ',obj.m_pepSeq, ' in ', obj.m_strSpecName, '!']);
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
end
