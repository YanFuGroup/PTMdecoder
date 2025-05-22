function review_score = reviewScoreing(obj)
% Give a weighted score to help the researcher finding good matches.

if isempty(obj.m_peptides)
    review_score = -1;
    return;
end

% if the score is out of range
if obj.m_main_score>obj.m_main_score_upper_bound || ...
        obj.m_main_score<obj.m_main_score_lower_bound
    error(['The main score of match for ', obj.m_spectrum, ' is out ' ...
        'of the range. Please check the score.']);
end

% preprocess the obj.m_spectrum
% filter the peaks whose intensity is less than 0.01 max intensity
% obj.m_spectrum.peaks = peakProcess(obj.m_spectrum.peaks, 0.01);

% score No.0, for trypsin, there cannot be modification on C-term amino acid
% the score is either 1 or -inf
score0_enzyme_specificity = 1;
enzyme = 'trypsin';
for idx_pep = 1:length(obj.m_peptides)
    if isequal(enzyme,'trypsin')
        idx_find = find(obj.m_peptides(idx_pep).mod_pos == length(obj.m_peptides(idx_pep).seq), 1);
        if ~isempty(idx_find) && ...
                obj.m_peptides(idx_pep).mod_mass(idx_find) < 14.0266
            % If a modification occurs in C-term, and the mass is greater
            % then 14.0266 (average mass of methylation, slightly
            % greater than the monoisotopic mass of methylation, which
            % is 14.015650), the PSM is absolutely wrong. 
            score0_enzyme_specificity = -inf;
        end
    end
end


% score No.1, normalized main score, is in [0, 1]
score1_norm_main_score = (obj.m_main_score-obj.m_main_score_lower_bound)/...
    (obj.m_main_score_upper_bound-obj.m_main_score_lower_bound);
if ~obj.m_is_high_score_better
    score1_norm_main_score = 1-score1_norm_main_score;
end


% score No.2, consecutive match score, is in [0, 1]
%   report the average consecutive match score when there are more than one
%   peptide matches the spectrum
score2_consec_score = 0;
score3_localization_score = 0;
% score No.3, localization score, is in [-0.75, 0.25]
for idx_pep = 1:length(obj.m_peptides)
    % the fifth column of the obj.m_all_match_ions means peptide index
    match_ions = obj.m_all_match_ions(obj.m_all_match_ions(:,5)==idx_pep,[1,2]);
    
    % match number and gap number generate score2
    match_type_pos = unique(match_ions(:,[1,2]),'rows');
    gap_num = count_gap(match_type_pos, length(obj.m_peptides(idx_pep).seq));
    score2_consec_score = score2_consec_score + (size(match_type_pos,1)-gap_num*0.5)/...
        (2 * (length(obj.m_peptides(idx_pep).seq)-1));
    
    % match modifications localization peaks
    mod_poses = obj.m_peptides(idx_pep).mod_pos;
    mod_poses(mod_poses==0) = 1;
    mod_poses(mod_poses==length(obj.m_peptides(idx_pep).seq)+1) = ...
        length(obj.m_peptides(idx_pep).seq);
    num_need_loc_pos = 2*length(mod_poses); % number of position needed to be localized
    num_need_loc_pos = num_need_loc_pos-sum(mod_poses==1);
    num_need_loc_pos = num_need_loc_pos-...
        sum(mod_poses==length(obj.m_peptides(idx_pep).seq));
    num_known_loc_pos = 0; % number of position localized by peak
    for idx_pos = 1:length(mod_poses)
        % left position
        % b localization ion (2,pos-1), y localization ion (1,length-pos+1)
        if any((match_ions(:,1)==2 & match_ions(:,2)==mod_poses(idx_pos)-1)) || ...
                (any(match_ions(:,1)==1 & match_ions(:,2)==...
                length(obj.m_peptides(idx_pep).seq)-mod_poses(idx_pos)+1))
            num_known_loc_pos = num_known_loc_pos+1;
        end
        % right position
        % b localization ion (2,pos), y localization ion (1,length-pos)
        if any((match_ions(:,1)==2 & match_ions(:,2)==mod_poses(idx_pos))) || ...
                (any(match_ions(:,1)==1 & match_ions(:,2)==...
                length(obj.m_peptides(idx_pep).seq)-mod_poses(idx_pos)))
            num_known_loc_pos = num_known_loc_pos+1;
        end
    end
    score3_localization_score = score3_localization_score + ...
        num_known_loc_pos / num_need_loc_pos;
end
% average the consecutive score among peptides
score2_consec_score = score2_consec_score/length(obj.m_peptides);
score3_localization_score = score3_localization_score/length(obj.m_peptides);
score3_localization_score = (score3_localization_score-0.75)/0.25;


% score No.4, high m/z interval unmatching high peak score, is in (-inf,0]
score4_unmatching_score = 1;
is_expe_unmatch = ones(size(obj.m_spectrum.peaks,1),1);
is_expe_unmatch(obj.m_all_match_ions(:,4)) = 0;
% the max peak in an interval should be greater than 10% max intensity,
% otherwise, the interval is negligible
peak_negl_thres = 0.1*max(obj.m_spectrum.peaks(:,2));
% check 100Da intervals from precursor-50Da to 200~300 Da
ratio_left = 0.6; % neglect the peak whose intensity is under (this*max_interval)
right_bound = obj.m_spectrum.pre_mz-50;
while right_bound>300
    left_bound = right_bound-100;
    is_interv_peak = obj.m_spectrum.peaks(:,1) >= left_bound & ...
        obj.m_spectrum.peaks(:,1) < right_bound;
    interval_peaks = obj.m_spectrum.peaks(is_interv_peak, :);
    val_max_interval = max(interval_peaks(:, 2));
    % only consider the intervals whose max peak is high enough
    if val_max_interval>peak_negl_thres
        % find the relative high peaks
        is_high_enough = is_expe_unmatch(is_interv_peak) & ...
            interval_peaks(:,2) > ratio_left*val_max_interval & ...
            interval_peaks(:,2) > peak_negl_thres;
        score4_unmatching_score = score4_unmatching_score - ...
            sum((interval_peaks(is_high_enough,2)-ratio_left*val_max_interval)/...
            ((1-ratio_left)*val_max_interval));
    end
    right_bound = left_bound;
end
% check 100Da intervals from precursor+50Da to [max] Da
ratio_right = 0.3; % neglect the peak whose intensity is under (this*max_interval)
left_bound = obj.m_spectrum.pre_mz+50;
while left_bound>max(obj.m_spectrum.peaks(:,1))
    right_bound = left_bound+100;
    is_interv_peak = obj.m_spectrum.peaks(:,1) >= left_bound & ...
        obj.m_spectrum.peaks(:,1) < right_bound;
    interval_peaks = obj.m_spectrum.peaks(is_interv_peak, :);
    val_max_interval = max(interval_peaks(:, 2));
    % only consider the intervals whose max peak is high enough
    if val_max_interval>peak_negl_thres
        % find the relative high peaks
        is_high_enough = is_expe_unmatch(is_interv_peak) & ...
            interval_peaks(:,2) > ratio_right*val_max_interval & ...
            interval_peaks(:,2) > peak_negl_thres;
        score4_unmatching_score = score4_unmatching_score - ...
            sum((interval_peaks(is_high_enough,2)-ratio_right*val_max_interval)/...
            ((1-ratio_right)*val_max_interval));
    end
    left_bound = right_bound;
end

% integrate the scores
% review_score = obj.m_weight_factor*...
review_score = [score0_enzyme_specificity;
    score1_norm_main_score;
    score2_consec_score;
    score3_localization_score;
    score4_unmatching_score];
end




% function [peaks]=peakProcess(peaks,alpha)
% % Filter the peaks whose intensity is less than (alpha * max intensity)
% % input:
% %   peaks:  [m/z1, intens1; m/z2, intens2;...]
% %   alpha:  the filter factor
% % output:
% %   peaks:  peaks after preprocess
% 
% filter_thres = alpha*max(peaks(:,2));
% peaks = peaks(peaks(:,2)>filter_thres,:);
% end




function gap_num = count_gap(match_ions, pep_length)
% Count the gap on peptide sequence
% input:
%   match_ions: [match_type, match_pos]
%       match_type: the length is equal to the number of matched
%           peaks, 1 means b-ion, 2 means y-ion
%       match_pos:  same length as match_type, each element shows
%           the position of the ion on sequence
%   pep_length: the length of the peptide sequence
% output:
%   gap_num:    the number of gap (missing ion on consecutive sequence)

ion_types = [1, 2]; % only consider the b/y ion here
gap_num = 0;
for idx_type = 1:length(ion_types)
    match_pos_on_seq = match_ions(match_ions(:,1)==ion_types(idx_type),2);
    last_pos = 0; % last viewed position
    for idx_pos = 1:length(match_pos_on_seq)
        if match_pos_on_seq(idx_pos) ~= last_pos+1
            % count a gap when skip any amino acid
            gap_num = gap_num+1;
        end
        last_pos = match_pos_on_seq(idx_pos);
    end
    if match_pos_on_seq(idx_pos) ~= pep_length-1
        % count a gap if the last possible cleavage do not appear
        gap_num = gap_num+1;
    end
end
end


