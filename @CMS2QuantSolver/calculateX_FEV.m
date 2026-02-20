function [X,ionTypePosCharge,ionIntens]=calculateX_FEV(vNonRedunTheoryIonMz,Peaks,case_OLS_intens_weight)
% Build design matrix X for FEV model.
% Input:
%   vNonRedunTheoryIonMz (L x T double)
%       Non-redundant ion table, each row is a fragment ion:
%       [m/z, type (1 is b ion, 2 is y ion), ion number (position), charge, 
%       number of modifications, class index, whether an IMP can generate this ion]
%   Peaks (K x 3 double)
%       Matched peaks [ion_index, normalized_intensity, raw_intensity].
%   case_OLS_intens_weight (char/string)
%       Weighting mode for X construction.
%           'none': no weighting (default)
%           'multi_self': weight by the sum of matched ion intensities per ion (self)
%           'multi_average_all': weight by the average of matched ion intensities per ion type (average across all ions of the same type)
%           'multi_average_self': weight by the average of matched ion intensities per ion type (average across ions of the same type that are present in the spectrum)
%           'multi_average_log': weight by the average of log-transformed matched ion intensities per ion type (average across all ions of the same type)
%           'multi_average_sqrt': weight by the average of sqrt-transformed matched ion intensities per ion type (average across all ions of the same type)
% Output:
%   X (N x P double)
%       Linear-system matrix for quantification.
%   ionTypePosCharge (U x 3 double)
%       Ion tuple [type, position, charge] included in right-half ion groups.
%   ionIntens (U x 1 double)
%       Aggregated ion intensity per ionTypePosCharge.
if nargin < 3 || isempty(case_OLS_intens_weight)
	case_OLS_intens_weight = 'none';
end

matPeaksBelong=vNonRedunTheoryIonMz(:,7:end);
IonTypes=vNonRedunTheoryIonMz(Peaks(:,1),2:4);

del_rows=[];
for iBelong=1:size(matPeaksBelong,1)
	if ~ismember(vNonRedunTheoryIonMz(iBelong,2:4),IonTypes,'rows')
		del_rows=[del_rows,iBelong]; %#ok<AGROW>
	end
end
matPeaksIntens=zeros(size(vNonRedunTheoryIonMz,1),max(vNonRedunTheoryIonMz(:,6)));
for iPeak=1:size(Peaks,1)
	matPeaksIntens(Peaks(iPeak,1),vNonRedunTheoryIonMz(Peaks(iPeak,1),6))=-Peaks(iPeak,2);
end
delIonKind=sum(matPeaksIntens)==0;  % Delete unmatched ion classes (columns on the right half)
matPeaksIntens(:,delIonKind)=[];
X=[matPeaksBelong,matPeaksIntens];

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
	weight = accumarray(ic,Peaks(:,2));
	[~, idx] = ismember(vNonRedunTheoryIonMz(:, 2:4), IonTypes, 'rows');
	nums_iontypes = accumarray(idx(idx > 0), 1, [size(IonTypes, 1), 1]);
	assert(length(weight) == length(nums_iontypes), 'Lengths of weight and nums_iontypes are not equal');
	weight = weight ./ nums_iontypes;

	sqrtWeights = zeros(size(idx));
	sqrtWeights(idx > 0) = sqrt(weight(idx(idx > 0)));
	X = X .* sqrtWeights;
elseif isequal(case_OLS_intens_weight, 'multi_average_self')
    % Weight the equation using relative ion intensity, take direct average, using only observed ions
	[IonTypes,~,ic] = unique(IonTypes,'rows');
	weight = accumarray(ic,Peaks(:,2));
	nums_iontypes = accumarray(ic, 1);
	weight = weight ./ nums_iontypes;
	[~, idx] = ismember(vNonRedunTheoryIonMz(:, 2:4), IonTypes, 'rows');
	sqrtWeights = zeros(size(idx));
	sqrtWeights(idx > 0) = sqrt(weight(idx(idx > 0)));
	X = X .* sqrtWeights;
elseif isequal(case_OLS_intens_weight, 'multi_average_log')
    % Weight the equation using relative ion intensity, take logarithmic average
	[IonTypes,~,ic] = unique(IonTypes,'rows');
	post_intens = log(Peaks(:,3));
	post_intens = post_intens / max(post_intens);
	weight = accumarray(ic,post_intens);
	[~, idx] = ismember(vNonRedunTheoryIonMz(:, 2:4), IonTypes, 'rows');
	nums_iontypes = accumarray(idx(idx > 0), 1, [size(IonTypes, 1), 1]);
	assert(length(weight) == length(nums_iontypes), 'Lengths of weight and nums_iontypes are not equal');
	weight = weight ./ nums_iontypes;
	sqrtWeights = zeros(size(idx));
	sqrtWeights(idx > 0) = sqrt(weight(idx(idx > 0)));
	X = X .* sqrtWeights;
elseif isequal(case_OLS_intens_weight, 'multi_average_sqrt')
    % Weight the equation using relative ion intensity, take square root average
	[IonTypes,~,ic] = unique(IonTypes,'rows');
	post_intens = sqrt(Peaks(:,3));
	post_intens = post_intens / max(post_intens);
	weight = accumarray(ic,post_intens);
	[~, idx] = ismember(vNonRedunTheoryIonMz(:, 2:4), IonTypes, 'rows');
	nums_iontypes = accumarray(idx(idx > 0), 1, [size(IonTypes, 1), 1]);
	assert(length(weight) == length(nums_iontypes), 'Lengths of weight and nums_iontypes are not equal');
	weight = weight ./ nums_iontypes;
	sqrtWeights = zeros(size(idx));
	sqrtWeights(idx > 0) = sqrt(weight(idx(idx > 0)));
	X = X .* sqrtWeights;
end

% Delete theoretical ion types (rows) that were not matched
X(del_rows,:) = [];

% Export ion metadata and intensity aligned with columns of right-half ion groups.
ionTypePosCharge=vNonRedunTheoryIonMz(:,[6,2:4]);
ionTypePosCharge=unique(ionTypePosCharge,'rows');
ionTypePosCharge(delIonKind,:)=[];
ionTypePosCharge=ionTypePosCharge(:,2:end);
ionIntens=sum(X,1);% Sum of intensities for each non-redundant ion type
ionIntens(1:size(vNonRedunTheoryIonMz,2)-6)=[];
ionIntens=ionIntens';
ionIntens=-ionIntens;
end