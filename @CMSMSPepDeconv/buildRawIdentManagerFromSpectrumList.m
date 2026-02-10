function rawIdentManager = buildRawIdentManagerFromSpectrumList(obj, spectrum_list)
% Build a raw identification manager from a spectrum list
% Input:
%   obj (CMSMSPepDeconv)
%       Processor instance
%   spectrum_list (struct array)
%       spectrum entries with fields:
%       - dataset_name
%       - spectrum_name
%       - peptidoform_list_str
%       - peptidoform_list_abun
% Output:
%   rawIdentManager (CIMPRawIdentManager)
%       per-raw identification store manager

rawIdentManager = CIMPRawIdentManager();
for idx_spec = 1:length(spectrum_list)
    dataset_name = spectrum_list(idx_spec).dataset_name;
    spectrum_name = spectrum_list(idx_spec).spectrum_name;
    peptidoform_strs = spectrum_list(idx_spec).peptidoform_list_str;
    peptidoform_abuns = spectrum_list(idx_spec).peptidoform_list_abun;
    [isorts, c_ref_isointens, c_mz, cur_ch] = obj.getProfiles(dataset_name, spectrum_name);
    lfMasses = get_masses_IMPs(peptidoform_strs,[obj.m_fixedModNameMass;obj.m_variableModNameMass]);
    rawStore = rawIdentManager.getOrCreate(dataset_name);
    rawStore = rawStore.appendSpecQuant(isorts, c_ref_isointens, c_mz, cur_ch, peptidoform_strs, lfMasses, peptidoform_abuns);
    rawIdentManager.setStore(dataset_name, rawStore);
end
end


% Get the mass of each IMPs
function lfMasses = get_masses_IMPs(cstrIMP,modNameMass)
% Input:
%   cstrIMP (K x 1 cellstr/string)
%       IMP names
%   modNameMass (M x 3 cell)
%       modification names, specificities, and masses
% Output:
%   lfMasses (K x 1 double) Da
%       masses of each IMP
lfMasses = zeros(length(cstrIMP),1);
for idx_imp = 1:length(cstrIMP)
    % Split the sequence of peptide and the modification
    mod_seq = cstrIMP{idx_imp};
    reg_exp = '\{(.*?)\}';
    [mod_str, seq_str] = regexp(mod_seq,reg_exp,'tokens','split');
    % Join the strings of sequence, delete the first and last "_" and count
    %   the masses.
    seq_str = strjoin(seq_str,'');
    seq_str([1,end]) = [];
    lfMasses(idx_imp) = sum(CConstant.vAAmass(seq_str-'A'+1));
    % Add the masses of modifications
    for idx_mod = 1:length(mod_str)
        is_notfound = true;
        for idx_mlist = 1:size(modNameMass,1)
            if isequal(modNameMass{idx_mlist,1},mod_str{idx_mod}{1})
                is_notfound = false;
                lfMasses(idx_imp) = lfMasses(idx_imp) + modNameMass{idx_mlist,3};
                break;
            end
        end
        if is_notfound
            error(['Unexpected modification is found: "',mod_str{idx_mod}{1},'"!']);
        end
    end
    % Add the mass of water
    lfMasses(idx_imp) = lfMasses(idx_imp) + CConstant.hmass*2 + CConstant.omass;
end
end
