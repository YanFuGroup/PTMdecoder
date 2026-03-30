function [isorts, c_ref_isointens, cur_mz, cur_ch] = getProfiles(cMgfDatasetIO, cMs12DatasetIO, cMsFileMapper, ms1_tolerance, mgf_name, spectrum_name)
% Read MS1 profile around one MS2 spectrum precursor.
% Inputs:
%    cMgfDatasetIO (CMgfDatasetIO)
%        MGF dataset reader.
%    cMs12DatasetIO (CMS12DatasetIO)
%        MS1/MS2 dataset reader.
%    cMsFileMapper (CMsFileMapper)
%        MGF-to-MS1 mapping helper.
%    ms1_tolerance (struct)
%        MS1 tolerance config with fields value and isppm.
%    mgf_name (1 x 1 char/string)
%        Dataset file name in MGF.
%    spectrum_name (1 x 1 char/string)
%        Spectrum name in MGF.
% Outputs:
%    isorts (1 x 1 double)
%        MS1 retention time of matched precursor scan.
%    c_ref_isointens (1 x 1 double)
%        Reference isotope intensity near precursor m/z.
%    cur_mz (1 x 1 double)
%        Precursor m/z from MGF.
%    cur_ch (1 x 1 double/int)
%        Precursor charge from MGF.

MS2ScanI = CMS2SpecNameUtils.parseMS2ScanNumber(spectrum_name);
[~, cur_ch, cur_mz] = cMgfDatasetIO.read_oneSpec(mgf_name, spectrum_name);

mgf_stem = erase(mgf_name, '.mgf');
ms12_stem = cMsFileMapper.get_ms1_stem(mgf_stem);

MS2_index = cMs12DatasetIO.m_mapNameMS2Index(ms12_stem);
idx_cur_scan = MS2_index(:, 2) == MS2ScanI;
MS1Scan = MS2_index(idx_cur_scan, 1);
MS1_index = cMs12DatasetIO.m_mapNameMS1Index(ms12_stem);
MS1_peaks = cMs12DatasetIO.m_mapNameMS1Peaks(ms12_stem);
index_starts_MS1 = [1; MS1_index(1:size(MS1_index, 1), 3)];
ino = find(MS1_index(:, 1) == MS1Scan);
isorts = MS1_index(ino, 2);
IX = index_starts_MS1(ino):index_starts_MS1(ino + 1) - 1;
mz = MS1_peaks(IX, 1);
inten = MS1_peaks(IX, 2);
if ms1_tolerance.isppm
    ptol = ms1_tolerance.value * cur_mz * 1e-6;
else
    ptol = ms1_tolerance.value;
end

c_ptol = min([ptol, 0.3]);
c_ref_isointens = max(inten(abs(mz - cur_mz) < c_ptol));
if isempty(c_ref_isointens)
    c_ref_isointens = 0;
end
end