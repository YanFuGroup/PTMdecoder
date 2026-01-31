function plotMatch(obj)
% Plot the spectra and mark the matched peaks in the spectrum
% Input:
%   obj (CReviewPSM)
%       Review instance

% Get relative intensity of all peaks
experimental_peaks = obj.m_spectrum.peaks;
experimental_peaks(:,2) = experimental_peaks(:,2)/max(experimental_peaks(:,2));

% Separate the peak indexes into 3 group (unmatched, single matched and 
%   multiple matched), according to obj.m_all_match_ions:
%   [match_type, match_pos, charge, (*)expe_which, pep_which]
%   the fourth attribution
[matched_peaks, ~, matched_peaks_idxs] = unique(obj.m_all_match_ions(:,4),'rows');
count_idx = zeros(size(matched_peaks,1),1); % the number of appearance of each peak (row)
for idx = 1:length(matched_peaks_idxs)
    count_idx(matched_peaks_idxs(idx)) = count_idx(matched_peaks_idxs(idx))+1;
end
single_peak_indexes = matched_peaks(count_idx==1);
multi_peak_indexes = matched_peaks(count_idx>1);

% match precursor -> fragment ion
if obj.m_tolerance.is_ppm
    mz_tol = obj.m_spectrum.pre_mz * tolerance.value / 1e6;
else
    mz_tol = obj.m_spectrum.pre_mz;
end
precursor_index = find(abs(experimental_peaks(:,1)-obj.m_spectrum.pre_mz)<mz_tol);
if length(precursor_index) > 1
    % The highest peak as precursor peak
    [~,precursor_index] = max(experimental_peaks(precursor_index,2));
end

% Plot 
figure
set(gcf,'position',[50,50,900,600], 'color','white')
set(gca,'LooseInset',get(gca,'TightInset'));
hold on
ion_type = {'y','b'};

% Plot the precursor peak
X_precursor = experimental_peaks(precursor_index,1); % m/z of peaks
Y_precursor = experimental_peaks(precursor_index,2); % intensity of peaks
stem(X_precursor,Y_precursor,'|g','LineWidth',0.75); % stem plot
text(X_precursor,Y_precursor+0.01,['[M]',repmat('+',1,obj.m_spectrum.pre_charge)],'Rotation',90,'HorizontalAlignment','left')

% Plot the unmatched peaks
unmatched_peak_indexes = 1:size(experimental_peaks,1);
unmatched_peak_indexes([single_peak_indexes;multi_peak_indexes;precursor_index]) = [];
X_unmatched = experimental_peaks(unmatched_peak_indexes,1); % m/z of peaks
Y_unmatched = experimental_peaks(unmatched_peak_indexes,2); % intensity of peaks
stem(X_unmatched,Y_unmatched,'|','Color','#9E9E9E','LineWidth',0.75)            % stem plot


% Plot the single matched peaks

X_single = experimental_peaks(single_peak_indexes,1);
Y_single = experimental_peaks(single_peak_indexes,2);
stem(X_single,Y_single,'|r','LineWidth',1.5)
for idx = 1:size(single_peak_indexes,1)  % labels of single matched peaks
    selected_ions = obj.m_all_match_ions(...
        find(obj.m_all_match_ions(:,4)==single_peak_indexes(idx)),:);
    ion_label = [num2str(selected_ions(5)), ':', ion_type{selected_ions(1)},...
        num2str(selected_ions(2)), repmat('+',1,selected_ions(3))];
    text(X_single(idx),Y_single(idx)+0.01,ion_label,'Rotation',90,'HorizontalAlignment','left')
%     text(X_single(idx),Y_single(idx)+0.01,ion_label,'HorizontalAlignment','center')
end

% Plot the multi matched peaks
hold on
X_single = experimental_peaks(multi_peak_indexes,1);
Y_single = experimental_peaks(multi_peak_indexes,2);
stem(X_single,Y_single,'|b','LineWidth',1.5)
for idx = 1:size(multi_peak_indexes,1)  % labels of multi matched peaks
    selected_ions = obj.m_all_match_ions(...
        find(obj.m_all_match_ions(:,4)==multi_peak_indexes(idx)),:);
    ion_label = [];
    for idx_pep = 1:2
        ion_label = [ion_label, num2str(selected_ions(idx_pep,5)), ':', ...
            ion_type{selected_ions(idx_pep,1)}, ...
            num2str(selected_ions(idx_pep,2)), ...
            repmat('+',1,selected_ions(idx_pep,3)),';'];
    end
    if size(selected_ions,1) > 2
        ion_label = [ion_label,'...'];
    else
        ion_label(end) = [];
    end
    text(X_single(idx),Y_single(idx)+0.01,ion_label,'Rotation',90,'HorizontalAlignment','left')
end
xlabel('m/z')
ylabel('Relative intensity (%)')
% saveas(gcf,'test.svg') % save into vector graph
end