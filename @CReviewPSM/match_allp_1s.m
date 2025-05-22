function obj = match_allp_1s(obj)
% Match all theoretical peptide ions with peaks in the spectrum.
% (output:)
%   obj.all_match_ions: ->  [match_type, match_pos, charge, expe_which, pep_which]
%       match_type: the length is equal to the number of matched 
%           peaks, 1 means b-ion, 2 means y-ion
%       match_pos:  same length as match_type, each element shows 
%           the position of the ion on sequence
%       charge:     the charge of theoretical ion
%       expe_which: the index of matched experimental peak in 
%           experimental spectrum
%       pep_which:  the index of matched peptide in peptides list

% all of the matched ions
obj.m_all_match_ions = zeros(2000,5);
i_all_match = 1;
for idx_pep = 1:length(obj.m_peptides)
    % peptide(s)-spectrum match
    match_ions = obj.match_1p_1s(obj.m_peptides(idx_pep), obj.m_spectrum, ...
        obj.m_tolerance); % may add types in the future
    if i_all_match+size(match_ions,1)>size(obj.m_all_match_ions,1)
        obj.m_all_match_ions(size(obj.m_all_match_ions,1)+1,...
            size(obj.m_all_match_ions,1)+2000,:) = zeros(2000,5);
    end
    obj.m_all_match_ions(i_all_match:i_all_match+size(match_ions,1)-1,:) = ...
        [match_ions,repmat(idx_pep,size(match_ions,1),1)];
    i_all_match = i_all_match+size(match_ions,1);
end
obj.m_all_match_ions(i_all_match:end,:) = [];
end

