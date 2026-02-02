function idxs_res = cluster_by_mass(imp_masses, ms1_tolerance)
    % Cluster the IMPs according to their masses.
    % Input:
    %   imp_masses (1 x K double) Da
    %   ms1_tolerance (struct) fields: isppm, value
    % Output:
    %   idxs_res (1 x G cell) indices of each group

    [m_val, m_inx] = sort(imp_masses);

    idxs_res = {};
    currentCluster = [];

    for i = 1:length(m_val)-1
        currentNumber = m_val(i);
        nextNumber = m_val(i+1);
        currentCluster = [currentCluster, m_inx(i)]; %#ok<AGROW>

        if ms1_tolerance.isppm
            tol = (ms1_tolerance.value * currentNumber)/1e6;
        else
            tol = ms1_tolerance.value;
        end

        if abs(nextNumber - currentNumber) > tol
            idxs_res = [idxs_res, {currentCluster}]; %#ok<AGROW>
            currentCluster = [];
        end
    end

    currentCluster = [currentCluster, m_inx(end)];
    idxs_res = [idxs_res, {currentCluster}];
end