function [bSuccess,cstrIMP,abundance,ionTypePosCharge,ionIntens,frageff,warning_msg,is_X_full_column_rank]=runEach(obj)
% Main entry point for discrimination and quantification of an IMP in a single spectrum
% Output: 
%   bSuccess - whether successful
%   cstrIMP - modification combination, character string form
%   abundance - content of various IMPs
%   ionTypePosCharge - fragment ion information involved in decomposition [by type, position, charge], each column corresponds to a column on the right half of the X matrix
%   ionIntens - intensity of each ion in ionTypePosCharge
%   frageff - fragmentation efficiency, arranged in the order of m_ionTypePosCharge
%       ionIntens and frageff are used to organize the fragmentation efficiency of ions involved in the model
%   warning_msg - the warning message
%   is_X_full_column_rank - true if X is full column rank, false otherwise

ionTypePosCharge=[];
ionIntens=[];
frageff=[];
is_X_full_column_rank=false;

%% Read required information from the spectrum
%   expPeaks - each row is an experimental spectral peak, left is m/z, right is intensity; 
%   iCharge - charge; 
%   dPrecursorMZ - neutral precursor ion mass
[obj.m_expPeaks,obj.m_iCharge,dPrecursorMZ]=obj.m_cMgfDatasetIO.read_oneSpec( ...
    obj.m_strDatasetName,obj.m_strSpecName);
obj.m_dPrecursorMass=(dPrecursorMZ-CConstant.pmass)*obj.m_iCharge;

%% Generate various modification configurations
% Get fixed modification information on the sequence, empty if none, otherwise each row is a fixed modification,
%   row content is [position, name (in Carbamidomethyl form), mass]
fixedPosMod=obj.getFixedPosMod();
[ bSuccess,inxSites,massArrangement,warning_msg ] = obj.getMassArrangement(fixedPosMod);
if ~bSuccess
    cstrIMP=[];
    abundance=[];
    return;
else
    %% Match multiple possible situations
    [vNonRedunTheoryIonMz]=obj.getNonRedunIons(inxSites,massArrangement,fixedPosMod);

    % Get peaks from the spectral dataset
    [expPeaks,~,~] = obj.m_cMgfDatasetIO.read_oneSpec(obj.m_strDatasetName,obj.m_strSpecName);%Each row is an experimental peak, left is m/z, right is intensity
    % expPeaks = peakPreprocess(expPeaks,obj.m_alpha);
    matchedExpPeaks = obj.match(expPeaks,vNonRedunTheoryIonMz); % Matching peaks using theoretical m/z of all configurations, remove unmatched peaks
    % matchedExpPeaks = normalize(matchedExpPeaks);
    matchedExpPeaks = peakProcess(matchedExpPeaks,obj.m_alpha); % Preprocess matched peaks

    % Cannot discriminate the result when there is no non-redundant peaks
    %   and several candidates.
    if isempty(matchedExpPeaks)
        if size(massArrangement,1) ~= 1
            cstrIMP = [];
            abundance = [];
            bSuccess = false;
            warning_msg = ['There is no non-redundant peak for discriminating the' ...
                ' peptidoforms for ',obj.m_pepSeq, ' in ', obj.m_strSpecName, '!\n'];
            return;
        else % Only one possible peptidoform
            abundance = 1;
            [cstrIMP]=obj.massArraTostr(massArrangement,fixedPosMod,...
                obj.m_variableModNameMass,inxSites);
            return;
        end
    end

    [massArrangement, vNonRedunTheoryIonMz] = obj.delete_useless_peptidoforms( ...
        matchedExpPeaks, massArrangement, vNonRedunTheoryIonMz);
    if size(massArrangement,1)==1
        abundance = 1;
        [cstrIMP]=obj.massArraTostr(massArrangement,fixedPosMod,...
            obj.m_variableModNameMass,inxSites);
        return;
    end


    %% Core model
    if obj.m_method == 3
        penalty_factor = obj.calculatePenaltyFactor(vNonRedunTheoryIonMz,matchedExpPeaks,obj.m_lambda);
    end
    switch obj.m_model
        case 1 % Fragmentation efficiency variable
            [X,ionTypePosCharge,ionIntens]=obj.calculateX_FEV(vNonRedunTheoryIonMz,matchedExpPeaks);    % A matrix
            switch obj.m_method
                case 1 % OLS
                    abundance=obj.coreFEV_OLS(X,massArrangement);
                case 2 % lasso
                    [abundance,frageff]=obj.coreFEV_lasso(X,massArrangement,obj.m_lambda);
                case 3 % penalty methods
                    [abundance,frageff]=obj.coreFEV_penalty(X,massArrangement,penalty_factor);
            end
        case 2 %Fragmentation efficiency constant - Guan
            X=obj.calculateX_Guan(vNonRedunTheoryIonMz);
            Y=obj.calculateY_Guan(vNonRedunTheoryIonMz,matchedExpPeaks);
            switch obj.m_method
                case 1 % OLS
                    abundance=obj.coreGuan_OLS(X,Y);
                case 2 % lasso
                    abundance=obj.coreGuan_lasso(X,Y,obj.m_lambda);
                case 3 % penalty methods
                    abundance = obj.coreGuan_penalty(X, Y, penalty_factor);
            end
        case 3 % Fragmentation efficiency equal
            X=obj.calculateX_Guan(vNonRedunTheoryIonMz);
            Y=obj.calculateY_1(vNonRedunTheoryIonMz,matchedExpPeaks);
            switch obj.m_method
                case 1 % OLS
                    abundance=obj.coreGuan_OLS(X,Y);
                case 2 % lasso
                    abundance=obj.coreGuan_lasso(X,Y,obj.m_lambda);
            end
    end

    % Determine if X is a singular matrix
    if rank(X)~=size(X,2)
        is_X_full_column_rank=true;
    end

    % Remove values less than the 1e-2 threshold, considered non-existent if smaller.
    abundance(abundance<obj.m_resFilterThres*max(abundance))=0;
    abundance=abundance/(sum(abundance)+eps);
end

% Interpret modification combinations as various IMP strings
[cstrIMP]=obj.massArraTostr(massArrangement,fixedPosMod,...
    obj.m_variableModNameMass,inxSites);

end



function [RelativePeaks]=peakProcess(Peaks,alpha)
% "Normalization" plus "relative intensity threshold denoising"
% Input: 
%   Peaks - the experimental spectrum, each row is a peak, left is m/z, right is intensity
%   alpha - relative intensity threshold, the peak with an intensity less than 'alpha' (e.g. 1%) of the maximum peak is considered noise and removed
% Output: 
%   RelativePeaks - the experimental spectrum, each row is a peak, first column is m/z, second column is normalized intensity, third column is original intensity

% Peaks(:,2)=log(Peaks(:,2));%Take natural logarithm ln of intensity first
% Peaks(:,2)=log10(Peaks(:,2));%Take base-10 logarithm log10 of intensity first
% Peaks(:,2)=sqrt(Peaks(:,2));%Take square root of intensity first
RelativePeaks=[Peaks,Peaks(:,2)];%Append original intensity to the third column
dFactor=max(RelativePeaks(:,2));
RelativePeaks(:,2)=RelativePeaks(:,2)/dFactor;%Normalize using the maximum value
RelativePeaks(RelativePeaks(:,2)<alpha,:)=[];%Remove noise peaks using relative intensity threshold denoising method
end