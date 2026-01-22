% Deprecated: use load_MS1_MS2_mapping instead

function success = load_MS2_file(obj,ms2_fullfile)
% Save the ms2 file as two files: "Spectrum Index" and "Peak Index"
% "Spectrum Index" (*_MS2_index.mat) stores the index of ten items such as spectrum name, as follows:
%   MS1 scan number, retention time, MS2 scan number, precursor m/z, charge state, fragmentation type
%   index of the first peak of the next spectrum, baseline for filtering peaks, MS2 retention time (same as the second item above, obtained in the same way), IonInjectionTime
% "Peak Index" (_MS2_peaks.mat) stores all peaks, with m/z on the left and intensity on the right, arranged densely without separation.
% Input: Name of the MS2 file, including path.
% Output: Success status, 1 for success, 0 for failure

success = 0;

% check the MAT file, if it exists, it means it has been processed before, skip processing and return success
[datapath,dataname] = fileparts(ms2_fullfile);
MS2_scanfile = fullfile(datapath,[dataname,'_MS2scans.mat']);
MS2_peakfile = fullfile(datapath,[dataname,'_MS2peaks.mat']);
if 0~=exist(MS2_scanfile,'file') && 0~=exist(MS2_peakfile,'file')
    success = 1;
    return;
end

% check the TXT file
if 0==exist(ms2_fullfile,'file')
    disp([ms2_fullfile,': does not exist!']);
    return;
end

% open the file
fid=fopen(ms2_fullfile,'r');
if -1==fid
    disp(['can not open: ',ms2_fullfile]);
    return;
end

%% init
% initialize the MS2 info
pmass = 1.007276;
maxpeaknum = 1e4; % the init peak number of a MS2 spectrum
maxMS2num = 1.5e5;% initial MS2 scan number
totalpeaknum = 4e7;% the init total peak number on MS2 scans
% MS1 scan,MS2 scan,m/z,z,MH,MS2 rt,Fragtype,IonInjectionTime,ActivationType(1-CID;2-ETD;3-HCD),InstrumentType(1-IT;2-FT),MS2 peak num
MS2_index = zeros([maxMS2num,10]);% MS1 scan,MS1 rt,MS2 scan,m/z,z,Fragtype,MS2 peak num,baseline,MS2 rt,MS2 injection time
MS2_peaks = zeros([totalpeaknum,2]);% m/z and intensity on MS2
fno = 0;% real MS2 scan number
pkno = 0;% real total peak number
oldms2scan = 0;

MS1_scanfile = fullfile(datapath,[dataname,'_MS1scans.mat']);
load(MS1_scanfile,'MS1_index');
num_MS1 = size(MS1_index,1);

% get the keywords
keyword0 = 'H	DataType';
keyword1 = 'S';
% keyword2 = 'I	RetTime';
% keyword3 = 'I	PrecursorInt';
% keyword4 = 'I	IonInjectionTime';
% keyword5 = 'I	ActivationType';
% keyword6 = 'I	PrecursorFile';
% keyword7 = 'I	PrecursorScan';
% keyword8 = 'I	InstrumentType';
% keyword9 = 'I	ActivationCenter';
keyword2 = 'I	NumberOfPeaks';
keyword3 = 'I	RetTime';
keyword4 = 'I	IonInjectionTime';
keyword5 = 'I	ActivationType';
keyword6 = 'I	InstrumentType';
keyword7 = 'I	PrecursorScan';
keyword8 = 'I	ActivationCenter';
keyword9 = 'I	MonoiosotopicMz';
keyword10 = 'Z';
len0 = length(keyword0);
len1 = length(keyword1);
len2 = length(keyword2);%#ok
len3 = length(keyword3);
len4 = length(keyword4);
len5 = length(keyword5);
len6 = length(keyword6);
len7 = length(keyword7);
len8 = length(keyword8);
len9 = length(keyword9);%#ok
len10 = length(keyword10);

%% get the MS2 info and files.
% get the datatype
str=fgets(fid);
while feof(fid)==0 && 0==strcmp( str(1:len0),keyword0 )
    str=fgets(fid);
end
MS2_datamode = str(len0+2);

if 1==strcmp('P',MS2_datamode)% DataType can be Centroid starting with C, or Profile starting with P
    fprintf(1,'MS2 is profile mode, convert to centroid mode first!\n');
    fclose(fid);
    return;
end

% for progress
ct_prt = 0;
fprintf(1,'MS2 scans: ');

% start to process
str = fgets(fid);
while 0==feof(fid)
    if 1==strcmp( str(1:len1),keyword1 )
        % 1.get the MS2 info
        tmp_datam = zeros(1,10);
        
        % MS2 scan MS2 scan m/z
        tmp_datam(1:3) = str2num(str(len1+2:end));%#ok
        
        % check the ms2scan
        if 0==oldms2scan-tmp_datam(2)
            str=fgets(fid);
            while feof(fid)==0 && 0==strcmp( str(1:len1),keyword1 )
                str=fgets(fid);
            end
            continue;
        else
            oldms2scan = tmp_datam(2);
        end
        
        % progress
        fno = fno + 1;
        fprintf(repmat('\b',[1,ct_prt]));
        ct_prt = fprintf('%i',fno);
        
        % OrigiNumberOfPeaks
        str=fgets(fid);%#ok
        %     if 1==strcmp( str(1:len2),keyword2 )
        %         tmp_datam(7) = str2num(str(len2+2:end));%#ok
        %     end;
        
        % RetTime
        str=fgets(fid);
        if 1==strcmp( str(1:len3),keyword3 )
            tmp_datam(6) = str2num(str(len3+2:end));%#ok
            MS2_RT=tmp_datam(6);
        end
        
        % IonInjectionTime
        str=fgets(fid);
        if 1==strcmp( str(1:len4),keyword4 )
            tmp_datam(8) = str2num(str(len4+2:end));%#ok
            MS2_InjecTime=tmp_datam(8);
        end
        
        % ActivationType
        str=fgets(fid);
        if 1==strcmp( str(1:len5),keyword5 )
            tmp_str = str(len5+2:len5+4);
            if 1==strcmp( tmp_str,'CID' )
                tmp_datam(9) = 1;
            elseif 1==strcmp( tmp_str,'ETD' )
                tmp_datam(9) = 2;
            elseif 1==strcmp( tmp_str,'HCD' )
                tmp_datam(9) = 3;
            end
        end
        
        % InstrumentType
        str=fgets(fid);
        if 1==strcmp( str(1:len6),keyword6 )
            tmp_str = str(len6+2:len6+5);
            if 1==strcmp( tmp_str,'ITMS' )
                tmp_datam(10) = 1;
            elseif 1==strcmp( tmp_str,'FTMS' )
                tmp_datam(10) = 2;
            end
        end
        % Fragtype
        tmp_datam(7) = (tmp_datam(9)-1)*2 + tmp_datam(10);% MS2dirs = {'CIDIT','CIDFT','ETDIT','ETDFT','HCDIT','HCDFT'};
        
        % MS1 scan
        str=fgets(fid);
        if 1==strcmp( str(1:len7),keyword7 )
            tmp_datam(1) = str2num(str(len7+2:end));%#ok
            p = find(MS1_index(1:num_MS1,1)<=tmp_datam(2));% All indices in the spectrum count that are smaller than it
            if 1==isempty(p)
                fno = fno - 1;
                continue;
            end
            if 1==isempty(tmp_datam(1)) || 0==tmp_datam(1)
                tmp_datam(1) = MS1_index(p(end),1);% MS1 scan number
            end
            tmp_datam(6) = MS1_index(p(end),2);% retention time
        end
        
        % Using this as center of m/z, the closer an isotopic peak is to
        % this center, the more confident we choose that peak to be the
        % monoisotopic peak. When Z. C. took over this code, this
        % block is not commented.
        % ActivationCenter
        str=fgets(fid);
%         if 1==strcmp( str(1:len8),keyword8 )
%             cen_mz = str2num(str(len8+2:end));%#ok
%         end
        
        % MonoiosotopicMz
        str=fgets(fid);
        % This block was once commented before Z. C. took over this code
        if 1==strcmp( str(1:len9),keyword9 )
            cen_mz = str2num(str(len9+2:end));%#ok
        end
        
        % z; MH
        str=fgets(fid);
        if 1==strcmp( str(1:len10),keyword10 )
            tmp_datam(4:5) = str2num(str(len10+2:end));%#ok 4 is charge, 5 is mass
            acc_mz = obj.get_acc_mz(cen_mz,tmp_datam(3),tmp_datam(4));% 3 is precursor m/z of the spectrum
            acc_mh = acc_mz*tmp_datam(4)-(tmp_datam(4)-1)*pmass;
            tmp_datam(3) = acc_mz;% corrected m/z
            tmp_datam(5) = acc_mh;% corrected M+
        end
        str=fgets(fid);
        while 1==strcmp( str(1:len10),keyword10 )% check for more charge?
            tmp_datam(4:5) = [0 0];
            str=fgets(fid);
        end
        
        % 2.read the MS2 data
        mz = zeros([1,maxpeaknum]);
        inten = zeros([1,maxpeaknum]);
        pnum = 0;
        while feof(fid)==0 && 0==strcmp( str(1:len1),keyword1 )
            if ~('0'<=str(1) && str(1)<='9')
                str = fgets(fid);
                continue;
            end
            pnum = pnum + 1;
            tmp = sscanf(str,'%f',2);
            mz(pnum) = tmp(1);
            inten(pnum) = tmp(2);
            str = fgets(fid);
        end
        if 1==feof(fid)
            pnum = pnum + 1;
            tmp = sscanf(str,'%f',2);
            mz(pnum) = tmp(1);
            inten(pnum) = tmp(2);
        end
        IX = find(inten>0);% keep the part with intensity greater than 0
        mz = mz(IX);
        inten = inten(IX);
        
        % 3.judge whether to centroid MS2 scan
        
        % 4.save the MS2 info and peaks
        if 1==isempty(mz) || mz(end)-mz(1)<10
            fno = fno - 1;
            continue;
        end
        if length(inten)>100
            baseline = GetBaseline(inten);
            IX = find(inten>=baseline);
            mz = mz(IX);
            inten = inten(IX);
        else
            baseline = 0;
        end
        npk = length(mz);
        
        x = [1 6 2 3 4 7];
        % 1 is MS1 scan number, 6 is retention time, 2 is MS2 scan number, 3 is precursor m/z, 4 is charge state, 7 is fragmentation type
        new_datam = tmp_datam(x);
        % new_datam is the above information, npk is the current number of peaks, baseline is the baseline for filtering peaks,
        % MS2_RT is MS2 retention time (same as the second item above, obtained in the same way), MS2_InjecTime is IonInjectionTime
        MS2_index(fno,1:10) = [new_datam npk baseline MS2_RT MS2_InjecTime];
        MS2_peaks(pkno+1:pkno+npk,1) = mz;
        MS2_peaks(pkno+1:pkno+npk,2) = inten;
        pkno = pkno + npk;% record the total number of spectra
    else
        str = fgets(fid);
    end
end
fclose(fid);
fprintf(repmat('\b',[1,ct_prt]));
fprintf('%i',fno);
fprintf(1,'\n');

%% save the MS2 info
% filter the empty values, keep valid records
if fno<maxMS2num
    IX = 1:fno;
    MS2_index = MS2_index(IX,:);% keep valid MS1 spectra
end
tmp = MS2_index(1:fno,7);
MS2_index(1:fno,7) = cumsum(tmp) + 1;%record the cumulative number of spectra, indicating the index number of the next spectrum
% if MS2_index(fno,2)>1000% 20*60=1200, if the retention time is too long, convert from seconds to minutes
% Convert the second column of MS2_index from seconds to minutes due to the current quantification method.
MS2_index(1:fno,2) = MS2_index(1:fno,2)/60; % 
% end

if pkno<totalpeaknum
    IX = 1:pkno;
    MS2_peaks = MS2_peaks(IX,:);
end

% save the results
save(MS2_scanfile,'MS2_index');
save(MS2_peakfile,'MS2_peaks');

success = 1;
end



function baseline = GetBaseline(inten)
% Get the baseline intensity for discarding peaks
% Input: intensity vector
% Output: baseline
% Method: Take logarithm, divide into steps, create histogram, take the center of the highest bar, restore baseline

loginten = log10(inten);
t = min(loginten):0.08:max(loginten);
[n,xout] = hist(loginten,t);

[tmp,idx] = max(n);%#ok
baseline = 10^xout(idx);
end

