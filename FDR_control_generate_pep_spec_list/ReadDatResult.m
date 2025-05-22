function [result2] = ReadDatResult(pathin)
fidin = fopen(pathin);

if -1==fidin
    disp(['Failed to open file ' pathin]);
    return;
end

% alphabet = ['A' 'B' 'C' 'D' 'E' 'F' 'G' 'H' 'I' 'J' 'K' 'L' 'M' 'N' 'O' 'P' 'Q' 'R' 'S' 'T' 'U' 'V' 'W' 'X' 'Y' 'Z'];
% lettercode = ['A' 'R' 'N' 'D' 'C' 'E' 'Q' 'G' 'H' 'I' 'L' 'K' 'M' 'F' 'P' 'S' 'T' 'W' 'Y' 'V'];
% monoisotopic = [71.03711 156.10111 114.04293 115.02694 103.00919 129.04259 128.05858 57.02146 137.05891 113.08406 113.08406 128.09496 131.04049 147.06841 97.05276 87.03203 101.04768 186.07931 163.06333 99.06841];

bufnum = 2000;
result1(bufnum) = struct('DatasetName','','Index1',0,'Calc_neutral_pepmass',0,'massdiff',0,'num_match_ions',0,'peptide','',...
    'modification','','modificationlocation','','ionscore',0,'protein','');
numrlt = 0;
ct_prt = 0;
s = fgetl(fidin);
head1 = 'IT_MODS=';
% head2 = 'queries=';
% % head3 = 'q1_p1=';
% head4 = '_p1=0';% Number of missed cleavage sites
% head5 = '_p1=1';
% head6 = '_p1=2';
% head44 = '_p1=3';
% head55 = '_p1=4';
% head66 = '_p1=5';
cleavage_tol = 5; % tolerance of number of cleavage
cleavage_pattern = string([repmat('_p1=',cleavage_tol+1,1),num2str((0:cleavage_tol)')]);
head7 = 'q1_p1=';
head8 = 'Content-Type: application/x-Mascot; name="proteins"';
head9 = 'MODS=';
head10 = 'FILE=';

% Read fixed modification information
while ~strncmp(s,head9,length(head9))
    s = fgetl(fidin);
end
fixmod = s(6:end);
if ~isempty(fixmod)
    fixmods = cell(1,2);
    S = regexp(fixmod,',','split');
    for i = 1:length(S)
        modinfo = regexp(char(S(i)),' ','split');
        fixmods{i,1} = char(S(i));
        modposition = char(modinfo(2));
        fixmods{i,2} = modposition(2:end-1);
    end
else
    fixmods = '';
end

% Read spectrum file name
while ~strncmp(s,head10,length(head10))
    s = fgetl(fidin);
end
S = regexp(s,'=','split');
S = regexp(S{end},'\','split');
datasetname = char(S(length(S)));


% Read variable modification information
while ~strncmp(s,head1,length(head1))
    s = fgetl(fidin);
end
itmods = regexp(s(9:end),',','split');

% while ~strncmp(s,head2,length(head2))
%     s = fgetl(fidin);
% end
% spectrumnum = s(9:end);

% find the first query "q1_p1"
while ~contains(s,head7)
    s = fgetl(fidin);
end
% tic;
% read all p1=...(specified num) before finding the sign of the end
% " Content-Type: application/x-Mascot; name="proteins" "
while ~strncmp(s,head8,length(head8))
    
    if ~contains(s,cleavage_pattern)
        s = fgetl(fidin);
        continue
    end
    
    % for example, 
    % q6176_p1=1,604.354416,
    % -0.000182,2,KSTGGK,8,
    % 02000000,11.74,
    % 0002000000000000000,0,0;"P68433":0:10:15:2,"P84228":0:10:15:2
    % meanings are:
    % q[query number]_p[index number]=[cleavage number],[neutral pep mass],
    % [mass shift],[num_matched_ions],[sequence],[most inten peak num used],
    % [modification on sequence],[ionscore],
    % [useless for now]
    S = regexp(s,';','split');
    info = char(S(1));
    S1 = regexp(info,',','split');
    
    index = regexp(char(S1(1)),'_','split');
    index = char(index(1));
    index = str2double(index(2:end));
    
    calcneutralpepmass = str2double(char(S1(2)));
    massdiff = str2double(char(S1(3)));
    nummatchions = str2double(char(S1(4)));
    peptide = char(S1(5));
    
    % variable modification
    mod = char(S1(7));
    if isempty(strrep(mod,'0',''))
        % has no variable modification
        modification = '-';
        modificationlocation = '-';
    else
        % record the modifications on each position of the sequence
        mod_idx_seq = mod-'0';
        mod_idx_seq(mod_idx_seq>10) = mod_idx_seq(mod_idx_seq>10)-7; % 'A'->10
        modification = char(join(itmods(mod_idx_seq(mod_idx_seq>0)),','));
        modificationlocation = char(join(string(find(mod_idx_seq)-1),','));
    end
    
    % Add fixed modifications
    % TODO: The program will go wrong when a terminal modificatino is set
    %   to be a fixed modification. Since here it is used to generate the
    %   pep-spec list, it seems no need to optimize it right now.
    if ~isempty(fixmods)
        fixPos = [];
        fixMod = [];
        Size = size(fixmods);
        for i = 1:Size(1)
            for j = 1:length(peptide)
                if peptide(j) == char(fixmods{i,2})
%                     if j == 0
%                         if strfind(fixmods{i,1},'N-term')
%                             fixPos = [fixPos,',',0];
%                             fixMod = [fixMod,',',char(fixmods{i,1})];
%                         else
%                             fixPos = [fixPos,',',1];
%                             fixMod = [fixMod,',',char(fixmods{i,1})];
%                         end
%                     else
                    fixPos = [fixPos,',',num2str(j)];
                    fixMod = [fixMod,',',char(fixmods{i,1})];
                    %                     end
                end
            end
            if contains(char(fixmods{i,1}),'N-term')
                fixPos = [fixPos,',','0'];
                fixMod = [fixMod,',',char(fixmods{i,1})];
            end
        end
        if ~isempty(fixPos)
            fixPos = fixPos(2:end);
            fixMod = fixMod(2:end);
            Pos = regexp(fixPos,',','split');
            Mod = regexp(fixMod,',','split');
        else
            Pos = '';
            Mod = '';
        end
        finalpos = '';
        finalmod = '';
        varmod = modification;
        varpos = modificationlocation;
        if modification == '-'
            varmod = '';
            varpos = '';
        end
        for num = 1:length(Pos)
            if ~contains(varpos,char(Pos(num)))
                finalpos = [finalpos,',',char(Pos(num))];
                finalmod = [finalmod,',',char(Mod(num))];
            end
        end
        if ~isempty(varpos) && ~isempty(finalpos)
            finalpos = [finalpos,',',varpos];
            finalmod = [finalmod,',',varmod];
            finalpos = finalpos(2:end);
            finalmod = finalmod(2:end);
            modification = finalmod;
            modificationlocation = finalpos;
        else
            if ~isempty(finalpos) && isempty(varpos)
                finalmod = finalmod(2:end);
                finalpos = finalpos(2:end);
                modification = finalmod;
                modificationlocation = finalpos;
            end
        end
    end
    
    ionscore = str2double(char(S1(8)));
    
    
    protein1 = char(S(2));
    S2 = regexp(protein1,'"','split');
    protein = '';
    for i = 1:(length(S2)-1)/2
        S = regexp(char(S2(2*i)),' ','split');
        protein = [protein ,',',char(S(1))];%#ok
    end
    protein = protein(2:end);
    
    s = fgetl(fidin);
    numrlt = numrlt+1;
    
    if numrlt>bufnum
        bufnum = bufnum+2000;
        result1(bufnum) = struct('DatasetName','','Index1',0,'Calc_neutral_pepmass',0,'massdiff',0,'num_match_ions',0,'peptide','',...
            'modification','','modificationlocation','','ionscore',0,'protein','');
    end
    
    result1(numrlt).DatasetName = datasetname;
    result1(numrlt).Index1 = index;
    result1(numrlt).Calc_neutral_pepmass = calcneutralpepmass;
    result1(numrlt).massdiff = massdiff;
    result1(numrlt).num_match_ions = nummatchions;
    result1(numrlt).peptide = peptide;
    result1(numrlt).modification = modification;
    result1(numrlt).modificationlocation = modificationlocation;
    result1(numrlt).ionscore = ionscore;
    result1(numrlt).protein = protein;
    
    if rem(numrlt,200)==0
        for j=1:ct_prt
            fprintf('\b');
        end
        ct_prt = fprintf('%i..',numrlt);
    end
    
    
end
% toc;
result1(numrlt+1:end)=[];

% tic;
result2(length(result1)) = struct('Site','','DatasetName','','Scan','','Spectrum','','Charge',0,...
    'Calc_neutral_pepmass',0,'precursor_neutral_mass',0,'massdiff',0,'num_match_ions',0,'peptide','',...
    'protein','','modification','','modificationlocation','','Score',0);

ct_prt = 0;


for i = 1:length(result1)
    head = ['Content-Type: application/x-Mascot; name="query',num2str(result1(i).Index1),'"'];
    %     tic;
    while ~contains(s,head)
        s = fgetl(fidin);
    end
    %     toc;
    fgetl(fidin);
    s = fgetl(fidin);
    spectrum = s(7:end);
    
    % Translate URL encoded characters to normal ASCII characters
    spectrum = strrep(spectrum,'%2e','.');
    spectrum = strrep(spectrum,'%20',' ');
    spectrum = strrep(spectrum,'%3a',':');
    spectrum = strrep(spectrum,'%3d','=');
    spectrum = strrep(spectrum,'%22','"');
    spectrum = strrep(spectrum,'%2c',',');
    spectrum = strrep(spectrum,'%2d','-');
    
    
    while ~contains(s,'charge=')
        s = fgetl(fidin);
    end
    charge = s(8:end-1);
    
    
    %scan
    if ~contains(spectrum,'.')
        scanpos = strfind(spectrum,':');
        scan = spectrum(scanpos+1:end);
    else
        scanpos = strfind(spectrum,'.');
        scan = spectrum(scanpos(1)+1:scanpos(2)-1);
    end
    
    result2(i).Site = pathin;
    result2(i).DatasetName = result1(1).DatasetName;
    result2(i).Scan = scan;
    result2(i).Spectrum = spectrum;
    %result2(i).Index = index;
    result2(i).Charge = charge;
    %result2(i).Index1 = result1(i).Index1;
    result2(i).Calc_neutral_pepmass = result1(i).Calc_neutral_pepmass;
    result2(i).precursor_neutral_mass = result2(i).Calc_neutral_pepmass + result1(i).massdiff;
    result2(i).massdiff = result1(i).massdiff;
    result2(i).num_match_ions = result1(i).num_match_ions;
    result2(i).peptide = result1(i).peptide;
    result2(i).protein = result1(i).protein;
    result2(i).modification = result1(i).modification;
    result2(i).modificationlocation = result1(i).modificationlocation;
    %     result2(i).ionscore = result1(i).ionscore;
    %
    result2(i).Score = result1(i).ionscore;
    
    
    
    if rem(i,200)==0
        for j=1:ct_prt
            fprintf('\b');
        end
        ct_prt = fprintf('%i..',i);
    end
end
% toc;



fclose(fidin);


end


