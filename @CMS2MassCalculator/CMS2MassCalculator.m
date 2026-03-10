classdef CMS2MassCalculator
    % Static mass and theoretical ion calculations for MS2

    methods (Static)
        fixedPosMod = getFixedPosMod(ctx)

        m_vTheoryMass = getNeutralPeptideTheoryMass(ctx, fixedPosMod)

        theoryMz = calculateIonMz(ctx, fixedPosMod)

        [bSuccess, inxSites, massArrangement] = getMassArrangement(ctx, fixedPosMod)

        vNonRedunTheoryIonMz = getNonRedunIons(ctx, inxSites, massArrangement, fixedPosMod)
    end
end
