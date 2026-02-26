function pep_quant = readFdrPeptides(obj, input_file_path, ms12DatasetIO, peptide_list, pep_quant)
% Parse filtered PSM results and build raw managers for target peptides
% Input:
%   obj (CMSMSPepDeconv)
%       Processor instance
%       File identifier of the filtered result file
%   input_file_path (1 x 1 char/string)
%       Path to the filtered result file
%   ms12DatasetIO (object)
%       MS1/MS2 dataset IO instance
%   peptide_list (1 x K cell)
%       Target peptide sequences
%   pep_quant (1 x K cell)
%       Per-peptide raw identification store managers
% Output:
%   pep_quant (1 x K cell)
%       Updated raw identification managers

if nargin < 3 || isempty(ms12DatasetIO)
    [obj, ~] = obj.ensureMs12DatasetIO();
    ms12DatasetIO = obj.m_cMs12DatasetIO;
end

FDRfilteredResults = CFdrFilteredResultIO.read(input_file_path);
entries = FDRfilteredResults.entries;

deps = struct( ...
    'ms12DatasetIO', ms12DatasetIO, ...
    'ms1_tolerance', obj.m_ms1_tolerance);

pep_quant = CPeptideRawIdentAssembler.buildFromFdrEntries(entries, peptide_list, pep_quant, deps);
end
