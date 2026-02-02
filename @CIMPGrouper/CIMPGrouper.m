classdef CIMPGrouper
    % Utilities for grouping IMPs

    methods (Static)
        % Cluster the IMPs according to their masses.
        idxs_res = cluster_by_mass(imp_masses, ms1_tolerance)
    end
end
