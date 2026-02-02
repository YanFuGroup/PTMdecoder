classdef CIMPGroup
    % Data model for a grouped IMP set for a specific raw and charge

    properties
        rawName
        impNames
        impMass
        ratio
        rts
        intensity
        charge
        lowMzBound
        highMzBound
        selectedCharge
        chargeGroupIdxs
        impRtRanges
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
