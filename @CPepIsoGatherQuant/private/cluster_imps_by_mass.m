function idxs_res = cluster_imps_by_mass(imp_masses, ms1_tolerance)
% Cluster the IMPs according to their masses.
% Input:
%   imp_masses (1 x K double) Da
%       Masses of IMPs
%   ms1_tolerance (struct)
%       Tolerance of MS1 (fields: isppm, value)
% Output:
%   idxs_res (1 x G cell)
%       Indices of each group, in cell form

[m_val, m_inx] = sort(imp_masses);

% Initialize variables
idxs_res = {};
currentCluster = [];

% Iterate through the numbers
for i = 1:length(m_val)-1
    currentNumber = m_val(i);
    nextNumber = m_val(i+1);
    % Add the current number to the current cluster
    currentCluster = [currentCluster, m_inx(i)]; %#ok<AGROW>

    % Mass tolerance
    if ms1_tolerance.isppm
        tol = (ms1_tolerance.value * currentNumber)/1e6;
    else
        tol = ms1_tolerance.value;
    end
    
    % Check if the next mass is beyond the threshold
    if abs(nextNumber - currentNumber) > tol
        % Add the current cluster to the list of clusters
        idxs_res = [idxs_res, {currentCluster}]; %#ok<AGROW> 
        
        % Reset the current cluster
        currentCluster = [];
    end
end

% Add the last number to the current cluster
currentCluster = [currentCluster, m_inx(end)];

% Add the last cluster to the list of clusters
idxs_res = [idxs_res, {currentCluster}];
end
