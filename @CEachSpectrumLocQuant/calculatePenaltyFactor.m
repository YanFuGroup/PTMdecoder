function penalty_factor = calculatePenaltyFactor(~, NonRedunTheoryIonMz, matchedExpPeaks, lambda, case_penalty_intens, grid_penalty_intens)
% Calculate the penalty factor for the core model according to the given
%   theoretical ions and the matched experimental peaks.
% Input:
%   vNonRedunTheoryIonMz (L x M double):
%       Site-discrimining ions, each row is a fragment ion:
%       [m/z, type (1 is b ion, 2 is y ion), ion number (position), charge, 
%       number of modifications, class index, whether an IMP can generate this ion]
%   matchedExpPeaks (K x 3 double):
%       List of matched experimental peaks. 
%       The first column is the matched vNonRedunTheoryIonMz number, the second column is the normalized intensity, and the third column is the original intensity.
%   lambda (1 x 1 double):
%       Weight for the penalty factor.
%   case_penalty_intens (1 x 1 char/string):
%       Mode for scoring/penalty (e.g., 'intens_sum', 'log_intens_sum', 'sqrt_intens_sum', 'exp_intens_sum', 'log_intens_ksdp', 'sqrt_intens_ksdp', 'exp_intens_ksdp', 'sqrt_intens_ksdp_sqrt', 'intens_hyperscore').
%   grid_penalty_intens (1 x 1 char/string):
%       Mode for intensity aggregation inside KSDP (e.g., 'intens_sum', 'log_intens_sum', 'sqrt_intens_sum', 'exp_intens_sum').
% Output:
%   penalty_factor (P x 1 double):
%       Penalty factor for the core model.

if nargin < 5 || isempty(case_penalty_intens)
    case_penalty_intens = 'intens_sum';
end
if nargin < 6 || isempty(grid_penalty_intens)
    grid_penalty_intens = 'intens_sum';
end

% Currently, only the simplest matching scoring function is used, which is SDP. In this example, the intensities of all fragment ions corresponding to the IMP are directly summed up
%   as the score of the IMP, and then the reciprocal is taken as the penalty factor.
%   More complex scoring functions can also be used, such as various scoring functions currently used by search engines, such as XCorr, hyperscore, etc.

% Get the number of IMPs
num_nrti = size(NonRedunTheoryIonMz, 2) - 6;

% Initialize the penalty factor
penalty_factor = zeros(num_nrti, 1);

% Iterate over each column representing whether the IMP produces the corresponding theoretical ion
for idx_form = 1:num_nrti
    ionIndices = find(NonRedunTheoryIonMz(:, idx_form+6) == 1);
    if isequal(case_penalty_intens, 'intens_sum') || ...
            isequal(case_penalty_intens, 'log_intens_sum') || ...
            isequal(case_penalty_intens, 'sqrt_intens_sum') || ...
            isequal(case_penalty_intens, 'exp_intens_sum')
        % Calculate the score for the current column (IMP)
        form_score = sdp_score(ionIndices, matchedExpPeaks, case_penalty_intens);
    elseif isequal(case_penalty_intens, 'log_intens_ksdp') || ...
            isequal(case_penalty_intens, 'sqrt_intens_ksdp') || ...
            isequal(case_penalty_intens, 'exp_intens_ksdp') || ...
            isequal(case_penalty_intens, 'sqrt_intens_ksdp_sqrt')
        [ion_tag, matched_peaks_intensity] = get_ion_tag_intensity(NonRedunTheoryIonMz, matchedExpPeaks, ionIndices);
        form_score = ksdp_score(ion_tag, matched_peaks_intensity, grid_penalty_intens);
        if isequal(case_penalty_intens, 'sqrt_intens_ksdp_sqrt')
            form_score = sqrt(form_score);
        end
    elseif isequal(case_penalty_intens, 'intens_hyperscore')
        form_score = hyperscore_score(NonRedunTheoryIonMz, matchedExpPeaks, ionIndices);
    else
        error('Unknown case_penalty_intens.');
    end

    % Store the result in a vector
    penalty_factor(idx_form) = 1 / form_score;
end

penalty_factor = penalty_factor ./ min(penalty_factor);
penalty_factor = penalty_factor * lambda;

end



function form_score = sdp_score(ionIndices, matchedExpPeaks, case_penalty_intens)
% Calculate the score for the given ionIndices according to the matched
%   experimental peaks.
% Input:
%   ionIndices (Q x 1 double/int):
%       Indices of the ions to be calculated.
%   matchedExpPeaks (K x 3 double):
%       List of matched experimental peaks.
%   case_penalty_intens (1 x 1 char/string):
%       Mode for intensity aggregation.
% Output:
%   form_score (1 x 1 double):
%       Score for the given ions.

% Initialize the sum of the current column
form_score = 0;

for i_ion = 1:length(ionIndices)    % row number of NonRedunTheoryIonMz
    idx_matchedExpPeaks = find(matchedExpPeaks(:, 1) == ionIndices(i_ion));

    if isempty(idx_matchedExpPeaks)
        continue;
    end
    
    % Sum these peak intensities
    if isequal(case_penalty_intens, 'log_intens_sum')
        peakValue = matchedExpPeaks(idx_matchedExpPeaks, 3);
        form_score = form_score + log(peakValue);
    elseif isequal(case_penalty_intens, 'sqrt_intens_sum')
        peakValue = matchedExpPeaks(idx_matchedExpPeaks, 3);
        form_score = form_score + sqrt(peakValue);
    elseif isequal(case_penalty_intens, 'exp_intens_sum')
        peakValue = matchedExpPeaks(idx_matchedExpPeaks, 2);
        form_score = form_score + exp(peakValue);
    else
        peakValue = matchedExpPeaks(idx_matchedExpPeaks, 2);
        form_score = form_score + sum(peakValue);
    end
end
end



function [ion_tag, matched_peaks_intensity] = get_ion_tag_intensity(NonRedunTheoryIonMz, matchedExpPeaks, ionIndices)
% Get the ion tag and matched peaks intensity for the given ionIndices.
% Input:
%   NonRedunTheoryIonMz (L x M double):
%       Theoretical non-redundant ions.
%   matchedExpPeaks (K x 3 double):
%       List of matched experimental peaks.
%   ionIndices (Q x 1 double/int):
%       Indices of the ions to be calculated.
% Output:
%   ion_tag (T x P double):
%       Ion tag for the given ions. [by&charge type] * [position]. 
%       Matched ions are marked as 1, otherwise 0.
%   matched_peaks_intensity (K x 1 double):
%       Intensity of the matched peaks for the given ions.

% Get unique values from NonRedunTheoryIonMz
unique_type = unique(NonRedunTheoryIonMz(:, [2, 4]), 'rows');
unique_pos = unique(NonRedunTheoryIonMz(:, 3));

% Create a matrix with unique values
ion_tag = zeros(size(unique_type, 1), size(unique_pos, 1));

% Get the matched peaks intensity
matched_peaks_intensity = zeros(size(matchedExpPeaks, 1), 1);

% % Get the squared maximum intensity of the matched peaks
% sqrt_max_intens = max(sqrt(matchedExpPeaks(:, 3)));

% Fill the matrix with vNonRedunTheoryIonMz values
for i_peak = 1:size(matchedExpPeaks, 1)
    if ~ismember(matchedExpPeaks(i_peak, 1), ionIndices)
        continue;
    end
    idx_ion = matchedExpPeaks(i_peak, 1);
    idx_type = unique_type(:, 1) == NonRedunTheoryIonMz(idx_ion, 2) & ...
        unique_type(:, 2) == NonRedunTheoryIonMz(idx_ion, 4);
    idx_pos = unique_pos == NonRedunTheoryIonMz(idx_ion, 3);
    ion_tag(idx_type, idx_pos) = 1;
    % matched_peaks_intensity(idx_ion) = sqrt(matchedExpPeaks(i_peak, 3))/sqrt_max_intens;
    matched_peaks_intensity(i_peak) = matchedExpPeaks(i_peak, 2);
end
end



function form_score = ksdp_score(ion_tag, matched_peaks_intensity, grid_penalty_intens)
% Calculate the score for the given ionIndices according to the matched
%   experimental peaks.
% Input:
%   ion_tag (T x P double):
%       Ion tag for the given ions. [by&charge type] * [position]. 
%       Matched ions are marked as 1, otherwise 0.
%   matched_peaks_intensity (K x 1 double):
%       Intensity of the matched peaks for the given ions.
%   grid_penalty_intens (1 x 1 char/string):
%       Mode for intensity aggregation.
% Output:
%   form_score (1 x 1 double):
%       Score for the given ions.

l_win = 5;
gamma = 0.9;
alfa = 0.5;
l1 = floor( (l_win -1) / 2 );
l2 = ceil( (l_win -1) / 2 );

nions = numel(ion_tag);
len = size(ion_tag,2);
ion_tag = cumsum(ion_tag<=0,2);
ion_tag = [zeros(size(ion_tag,1),1) ion_tag];
kernel = 0;
for i=2:len+1
    % ion tag is a boolean matrix indicating the matched ions
    win = ion_tag(:,min(len+1,i+l2))-ion_tag(:,max(1,i-l1-1));
    kernel = kernel + sum(exp((-1)*gamma*win));
end

% Intensity sum in scoring function
if isequal(grid_penalty_intens, 'log_intens_sum')
    s_inten = sum(log(matched_peaks_intensity));
elseif isequal(grid_penalty_intens, 'sqrt_intens_sum')
    s_inten = sum(sqrt(matched_peaks_intensity));
elseif isequal(grid_penalty_intens, 'exp_intens_sum')
    s_inten = sum(exp(matched_peaks_intensity));
else
    s_inten = sum(matched_peaks_intensity);
end

form_score = s_inten * power(kernel,alfa) / nions;
end



function form_score = hyperscore_score(NonRedunTheoryIonMz, matchedExpPeaks, ionIndices)
% Calculate the score for the given peptidoforms according to the matched
%   experimental peaks.
% Input:
%   NonRedunTheoryIonMz (L x M double):
%       Theoretical non-redundant ions.
%   matchedExpPeaks (K x 3 double):
%       List of matched experimental peaks.
%   ionIndices (Q x 1 double/int):
%       Indices of the ions to be calculated.
%       According to the row numbers of NonRedunTheoryIonMz.
% Output:
%   form_score (1 x 1 double):
%       Score for the given peptidoforms.
% Attention:
%   The hyperscore is calculated mainly based on the b/y ions!

% hyperScore = log(Nb! * Ny! * \Sum(Inten_b) * \Sum(Inten_y))
%            = log(Nb!) + log(Ny!) + log(\Sum(Inten_b)) + log(\Sum(Inten_y)
%            Look up table   Look up table   Accumulate and take logarithm   Accumulate and take logarithm

form_score = 0;

% Limit the number of matched peaks, 1000 by default
MAX_MATCH_ION_NUM = 1000;
log_factorial = cumsum(log10(1:MAX_MATCH_ION_NUM));

% Get the numbers and intensities of matched b/y ions
Nb = 0;
Ny = 0;
sumInten_b = 0;
sumInten_y = 0;
for i_mp = 1:size(matchedExpPeaks, 1)
    if ~ismember(matchedExpPeaks(i_mp, 1), ionIndices)
        continue;
    end
    idx_ion = matchedExpPeaks(i_mp, 1);
    if NonRedunTheoryIonMz(idx_ion, 2) == 1
        Nb = Nb + 1;
        sumInten_b = sumInten_b + matchedExpPeaks(i_mp, 3);
    elseif NonRedunTheoryIonMz(idx_ion, 2) == 2
        Ny = Ny + 1;
        sumInten_y = sumInten_y + matchedExpPeaks(i_mp, 3);
    end
    if matchedExpPeaks(i_mp, 3) < 1
        error('The absolute intensity of the matched peak should be greater than 1.');
    end
    
    if Nb > MAX_MATCH_ION_NUM || Ny > MAX_MATCH_ION_NUM
        error('The number of matched peaks should be less than %d.', MAX_MATCH_ION_NUM);
    end

    if Nb > 0
        form_score = form_score + log_factorial(Nb) + log(sumInten_b);
    end
    if Ny > 0
        form_score = form_score + log_factorial(Ny) + log(sumInten_y);
    end
end
end