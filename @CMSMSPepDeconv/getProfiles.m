function [isorts,c_ref_isointens,cur_mz,cur_ch] = getProfiles(obj,mgf_name,spectrum_name)
% Get retention time and intensity of several (four) isotope peaks according to the specified MS1 spectrum. All are 0 when isotopes are incomplete
% input:    strDatasetName
%               name of dataset, mgf
%           strSpecName
%               name of spectra, dta
% output:   isorts
%               retention time
%           c_ref_isointens
%               vector of isotopic intensities
%           mz
%               precursor m/z of this spectra
%           cur_ch
%               precursor charge of this spectra

% ATTENTION: The second digit separated by a dot is required to be the scan number of the spectrum. If this condition is not met, a serious error will occur!!!
spec_name = regexp(spectrum_name,'\.','split');
MS2ScanI = str2double(spec_name{2});
[~, cur_ch, cur_mz] = obj.m_cMgfDatasetIO.read_oneSpec(mgf_name,spectrum_name);
MS2_index = obj.m_cMs12DatasetIO.m_mapNameMS2Index(erase(mgf_name,'.mgf'));
idx_cur_scan = MS2_index(:,2)==MS2ScanI; % Find the corresponding scan
MS1Scan = MS2_index(idx_cur_scan,1); % Record the scan number of MS1
% MS1_index (scan, retention time, peak number, baseline, injection time)
% MS1_peaks (m/z, intensity)
MS1_index = obj.m_cMs12DatasetIO.m_mapNameMS1Index(erase(mgf_name,'.mgf'));
MS1_peaks = obj.m_cMs12DatasetIO.m_mapNameMS1Peaks(erase(mgf_name,'.mgf'));
index_starts_MS1 = [1;MS1_index(1:size(MS1_index,1),3)]; % Starting position of each spectrum
ino = find(MS1_index(:,1)==MS1Scan); % Find the position of the MS1 scan number
isorts = MS1_index(ino,2); % Get its retention time
% Get the index range of the spectral peak, then obtain the m/z and intensity of the spectral peak within the index range, and record them in mz and inten respectively
IX = index_starts_MS1(ino):index_starts_MS1(ino+1)-1;
mz = MS1_peaks(IX,1);
inten = MS1_peaks(IX,2);
% tolerance
if obj.m_ms1_tolerance.isppm
    ptol = obj.m_ms1_tolerance.value*cur_mz*1e-6;
else
    ptol = obj.m_ms1_tolerance.value;
end
c_ptol = min([ptol,0.3]); % The final tolerance
c_ref_isointens = max(inten(abs(mz-cur_mz)<c_ptol));
if isempty(c_ref_isointens)
    c_ref_isointens = 0;
end
end
