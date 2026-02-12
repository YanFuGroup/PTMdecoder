classdef CIMPGroup
    % Data model for a grouped IMP set for a specific raw and charge

    properties
        rawName         % Raw file name for this group
        impNames        % IMP names in this group
        impMass         % IMP masses (Da)
        ratio           % IMP ratios per observation
        rts             % Retention times per observation (minutes)
        intensity       % Intensities per observation
        charge          % Charge state per observation
        lowMzBound      % Lower m/z bound for XIC extraction
        highMzBound     % Upper m/z bound for XIC extraction
        selectedCharge  % Charge used for alignment/quant
        chargeGroupIdxs % Indices of observations for selectedCharge
        impRtRanges     % Optional RT ranges per IMP (requant)
    end

    methods
        function obj = CIMPGroup(rawName, impNames, impMass, ratio, rts, intensity, charge, ...
                lowMzBound, highMzBound, selectedCharge, chargeGroupIdxs, impRtRanges)
            if nargin == 0
                return;
            end
            obj.rawName = rawName;
            obj.impNames = impNames;
            obj.impMass = impMass;
            obj.ratio = ratio;
            obj.rts = rts;
            obj.intensity = intensity;
            obj.charge = charge;
            obj.lowMzBound = lowMzBound;
            obj.highMzBound = highMzBound;
            obj.selectedCharge = selectedCharge;
            obj.chargeGroupIdxs = chargeGroupIdxs;
            obj.impRtRanges = impRtRanges;
        end
    end
end
