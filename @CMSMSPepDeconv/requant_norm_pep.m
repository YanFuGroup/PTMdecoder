function obj = requant_norm_pep(obj)
% Re-quantify the normalization peptides using checked XIC peaks
% Input:
%   obj (CMSMSPepDeconv)
%       Processor instance
% Output:
%   obj (CMSMSPepDeconv)
%       Updated instance

%% Read the checked peptides and their XIC range
if isempty(obj.m_checked_peptides_res_path)
    checked_pep_path = fullfile(obj.m_outputDir, 'peptide4normalization_checked.txt');
else
    checked_pep_path = obj.m_checked_peptides_res_path;
end

report = CIMPQuantResultIO.read(checked_pep_path);
pep_rtrange_map = report.build_pep_rtrange_map();
pep_prot_map = report.build_pep_prot_map();
if isempty(pep_rtrange_map)
    error(['The checked normalization peptide result file "', checked_pep_path, '" is empty!']);
end

if ~isfolder(obj.m_outputDir)
    mkdir(obj.m_outputDir);
end

%% Requantify the normalization peptides
output_path = fullfile(obj.m_outputDir, 'peptide4normalization_requant.txt');
report_requant = CIMPQuantReport();

% Indexing the dataset IO
obj.m_cMs12DatasetIO = CMS12DatasetIO(obj.m_specPath,obj.m_ms1_tolerance);
obj.m_cMs12DatasetIO.SetMap();
obj.m_cMgfDatasetIO = CMgfDatasetIO;
obj.m_cMgfDatasetIO.Init(obj.m_specPath);
obj.m_cMgfDatasetIO.SetMap();
obj.m_cMgfDatasetIO.SetFidmap();

% Read and process
total_records = 0;
for idx_block = 1:numel(report.blocks)
    total_records = total_records + numel(report.blocks(idx_block).records);
end
if total_records == 0
    fprintf(['Warning: The file "', checked_pep_path, '" is empty\n']);
end
print_progress = CPrintProgress(max(total_records, 1));
pipeline = CIMPProcessingPipeline(obj.m_cMs12DatasetIO, obj.m_ms1_tolerance, obj.m_min_MSMS_num, obj.m_alpha, obj.m_resFilterThres);

fprintf('Re-quantifying at peptide level...')
rec_counter = 0;
for idx_block = 1:numel(report.blocks)
    prot_name = '';
    if ~isempty(report.blocks(idx_block).protein_name_pos)
        prot_name = report.blocks(idx_block).protein_name_pos{1,1};
    end
    records = report.blocks(idx_block).records;
    for idx_rec = 1:numel(records)
        rec_counter = rec_counter + 1;
        print_progress = print_progress.update_show(rec_counter);

        mgf_name = records(idx_rec).raw_name;
        current_charge = records(idx_rec).charge;
        current_peptide = records(idx_rec).imp_name;
        lfMass = get_mass_peptide(current_peptide);
        lfMz = (lfMass + current_charge * CConstant.pmass) / current_charge;
        current_key = CIMPQuantRecord.build_id(current_peptide, current_charge, mgf_name);
        if pep_prot_map.isKey(current_key)
            prot_name = pep_prot_map(current_key);
        end
        rawIdentManager = CIMPRawIdentManager();

        rt_peaks = records(idx_rec).rt_peaks;
        for idx_peak = 1:numel(rt_peaks)
            rt_median = (rt_peaks(idx_peak).rt_start + rt_peaks(idx_peak).rt_end) / 2;
            rawStore = rawIdentManager.getOrCreate(mgf_name);
            rawStore = rawStore.appendSpecQuant(rt_median, 1, lfMz, current_charge, {current_peptide}, lfMass, 1);
            rawIdentManager.setStore(mgf_name, rawStore);
        end

        block = pipeline.requantifyPeptideBlock({prot_name, -1}, rawIdentManager, pep_rtrange_map);
        report_requant = report_requant.append_block(block);
    end
end

print_progress.last_update();
fprintf('done.\n');

CIMPQuantResultIO.write(report_requant, output_path);
obj.m_cMgfDatasetIO.CloseAllFile();
end



%% Other functions

% Get the mass of each IMPs
function lfMass = get_mass_peptide(pep_seq)
% Input:
%   pep_seq (1 x 1 char/string)
%       peptide sequence
% Output:
%   lfMass (1 x 1 double) Da
%       mass of the peptide

% Add the mass of each amino acid
lfMass = sum(CConstant.vAAmass(pep_seq-'A'+1));

% Add the mass of water
lfMass = lfMass + CConstant.hmass*2 + CConstant.omass;
end
