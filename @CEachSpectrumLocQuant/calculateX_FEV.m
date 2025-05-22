function [X,ionTypePosCharge,ionIntens]=calculateX_FEV(~,vNonRedunTheoryIonMz,Peaks)
% Calculate the X matrix in $Y=X\alpha+\epsilon$
% Input: 
%   vNonRedunTheoryIonMz - Site-discrimining ions, each row is a fragment ion:
%       [m/z, type (1 is b ion, 2 is y ion), ion number (position), charge, 
%       number of modifications, class index, whether an IMP can generate this ion]
%   massArrangement - the various mass arrangements of modifications on the peptide
%   Peaks - the matched experimental spectrum peaks, the first column is the index matched to vNonRedunTheoryIonMz, the second column is the normalized intensity, the third column is the real intensity obtained experimentally
%   ms2_tolerance - the fragment ion matching error
% Output: 
%   X - the X matrix in $Y=X\alpha+\epsilon$
%   ionTypePosCharge - the ion information contained in the matrix, each row corresponds to a column on the right half of the X matrix
%   ionIntens - the matched intensity, organizes ion intensities to determine which ions to use and how, can be greater than 1
matPeaksBelong=vNonRedunTheoryIonMz(:,7:end); %Left half
IonTypes=vNonRedunTheoryIonMz(Peaks(:,1),2:4);%Matched ion classes

del_rows=[];  
for iBelong=1:size(matPeaksBelong,1)
    if ~ismember(vNonRedunTheoryIonMz(iBelong,2:4),IonTypes,'rows')
        del_rows=[del_rows,iBelong]; %#ok<AGROW> 
    end
end
matPeaksIntens=zeros(size(vNonRedunTheoryIonMz,1),max(vNonRedunTheoryIonMz(:,6)));  %Right half
for iPeak=1:size(Peaks,1)
    matPeaksIntens(Peaks(iPeak,1),vNonRedunTheoryIonMz(Peaks(iPeak,1),6))=-Peaks(iPeak,2);
end
delIonKind=sum(matPeaksIntens)==0;  % Delete unmatched ion classes (columns on the right half)
matPeaksIntens(:,delIonKind)=[];
X=[matPeaksBelong,matPeaksIntens];


global case_OLS_intens_weight
if isequal(case_OLS_intens_weight, 'multi_self')
    % Weight the equation using relative ion intensity
    weight = abs(sum(X(:,size(matPeaksBelong,2)+1:end),2));
    for i_w = 1:length(weight)
        if weight(i_w) ~= 0
            X(i_w,:) = X(i_w,:)*sqrt(weight(i_w));
        end
    end
elseif isequal(case_OLS_intens_weight, 'multi_average_all')
    % Weight the equation using the average of relative ion intensity
    [IonTypes,~,ic] = unique(IonTypes,'rows');
    weight = accumarray(ic,Peaks(:,2)); % Weight of each ion type
    % Get matching indices
    [~, idx] = ismember(vNonRedunTheoryIonMz(:, 2:4), IonTypes, 'rows');
    % Count occurrences of each IonType
    nums_iontypes = accumarray(idx(idx > 0), 1, [size(IonTypes, 1), 1]);% Number of each ion type
    assert(length(weight) == length(nums_iontypes), 'Lengths of weight and nums_iontypes are not equal');
    weight = weight ./ nums_iontypes; % Average
    
    sqrtWeights = zeros(size(idx));
    % Calculate the square root of weights for non-zero indices
    sqrtWeights(idx > 0) = sqrt(weight(idx(idx > 0)));
    % Update the X matrix
    X = X .* sqrtWeights;
elseif isequal(case_OLS_intens_weight, 'multi_average_self')
    % Weight the equation using relative ion intensity, take direct average, using only observed ions
    [IonTypes,~,ic] = unique(IonTypes,'rows');
    weight = accumarray(ic,Peaks(:,2)); % Weight of each ion type, sum of intensities first
    % Count occurrences of each IonType, considering only observed peaks
    nums_iontypes = accumarray(ic, 1);% Number of each ion type
    weight = weight ./ nums_iontypes; % Then calculate the average
    % Get matching indices
    [~, idx] = ismember(vNonRedunTheoryIonMz(:, 2:4), IonTypes, 'rows');
    % Create a sqrtWeights array of the same size as idx, initialized to zero
    sqrtWeights = zeros(size(idx));
    % Calculate the square root of weights for non-zero indices
    sqrtWeights(idx > 0) = sqrt(weight(idx(idx > 0)));
    % Update the X matrix
    X = X .* sqrtWeights;
elseif isequal(case_OLS_intens_weight, 'multi_average_log')
    % Weight the equation using relative ion intensity, take logarithmic average
    [IonTypes,~,ic] = unique(IonTypes,'rows');
    post_intens = log(Peaks(:,3));
    post_intens = post_intens / max(post_intens);
    weight = accumarray(ic,post_intens);        % Weight of each ion type, sum of intensities first
    % Get matching indices
    [~, idx] = ismember(vNonRedunTheoryIonMz(:, 2:4), IonTypes, 'rows');
    % Count occurrences of each IonType
    nums_iontypes = accumarray(idx(idx > 0), 1, [size(IonTypes, 1), 1]);% Number of each ion type
    assert(length(weight) == length(nums_iontypes), 'Lengths of weight and nums_iontypes are not equal');
    weight = weight ./ nums_iontypes; % Then calculate the average
    % Create a sqrtWeights array of the same size as idx, initialized to zero
    sqrtWeights = zeros(size(idx));
    % Calculate the square root of weights for non-zero indices
    sqrtWeights(idx > 0) = sqrt(weight(idx(idx > 0)));
    % Update the X matrix
    X = X .* sqrtWeights;
elseif isequal(case_OLS_intens_weight, 'multi_average_sqrt')
    % Weight the equation using relative ion intensity, take square root average
    [IonTypes,~,ic] = unique(IonTypes,'rows');
    post_intens = sqrt(Peaks(:,3));
    post_intens = post_intens / max(post_intens);
    weight = accumarray(ic,post_intens);        % Weight of each ion type, sum of intensities first
    % Get matching indices
    [~, idx] = ismember(vNonRedunTheoryIonMz(:, 2:4), IonTypes, 'rows');
    % Count occurrences of each IonType
    nums_iontypes = accumarray(idx(idx > 0), 1, [size(IonTypes, 1), 1]);% Number of each ion type
    assert(length(weight) == length(nums_iontypes), 'Lengths of weight and nums_iontypes are not equal');
    weight = weight ./ nums_iontypes; % Then calculate the average
    % Create a sqrtWeights array of the same size as idx, initialized to zero
    sqrtWeights = zeros(size(idx));
    % Calculate the square root of weights for non-zero indices
    sqrtWeights(idx > 0) = sqrt(weight(idx(idx > 0)));
    % Update the X matrix
    X = X .* sqrtWeights;
end

% Delete theoretical ion types (rows) that were not matched
X(del_rows,:) = [];

ionTypePosCharge=vNonRedunTheoryIonMz(:,[6,2:4]);% Record non-redundant ion types, the first column is the ion number, which must not be missing, otherwise the position correspondence will be incorrect
ionTypePosCharge=unique(ionTypePosCharge,'rows');
ionTypePosCharge(delIonKind,:)=[];
ionTypePosCharge=ionTypePosCharge(:,2:end);
ionIntens=sum(X,1);% Sum of intensities for each non-redundant ion type
ionIntens(1:size(vNonRedunTheoryIonMz,2)-6)=[];
ionIntens=ionIntens';
ionIntens=-ionIntens;
end