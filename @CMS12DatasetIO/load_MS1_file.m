function success = load_MS1_file(~,ms1_fullfile)
% Function: Save the ms1 file as two files: "Spectrum Index" and "Peak Index", and store them in the path of the ms1 file.
% "Spectrum Index" (_MS1peaks.mat) stores the index of five items such as spectrum name, as follows:
%   Spectrum starting position in peaks, retention time, number of peaks obtained after filtering low-abundance peaks,
%   baseline used for filtering peaks, IonInjectionTime in the MS1 file.
% "Peak Index" (_MS1scans.mat) stores all peaks, with m/z on the left and intensity on the right, arranged densely without separation.
% Input: ms1_fullfile is the name of the MS1 file, including path.
% Output: success indicates whether it is successful, 1 for success, 0 for failure

success = 0;

% check the MAT file, if it exists, it means it has been processed before, skip processing and return success
[datapath,dataname] = fileparts(ms1_fullfile);
MS1_scanfile = fullfile(datapath,[dataname,'_MS1scans.mat']);
MS1_peakfile = fullfile(datapath,[dataname,'_MS1peaks.mat']);
if 0~=exist(MS1_scanfile,'file') && 0~=exist(MS1_peakfile,'file')
    %{%for old data do not have 'MS1Type'
    load(MS1_scanfile);
    if 0==exist('MS1Type')%#ok
        MS1Type = 'FTMS';
        save(MS1_scanfile,'MS1_index','MS1Type');
    end
    %}
    success = 1;
    return;
end

% check the TXT file
if 0==exist(ms1_fullfile,'file')
    disp([ms1_fullfile,': does not exist!']);
    return;
end

% open the TXT file
fid = fopen(ms1_fullfile,'r');
if -1==fid
    disp([ms1_fullfile,': can not open!']);
    return;
end

%% init
% initialize the MS1 info
maxpeaknum = 1e4;% the max peak number on a MS1 scan
maxMS1num = 1.5e5;% initial MS1 scan number
totalpeaknum = 4e7;% the init total peak number on MS1 scans
MS1_index = zeros([maxMS1num,5]);% MS1 scan, MS1 rt, MS1 peak num, baseline, injection time
MS1_peaks = zeros([totalpeaknum,2]);% m/z and intensity on MS1
fno = 0;% real MS1 scan number
pkno = 0;% real total peak number

% get the keywords
keyword0 = 'H	DataType';
keyword1 = 'S';% the keyword to start record
keyword2 = 'I	RetTime';% the keyword to start record
keyword3 = 'I	InstrumentType';
keyword4='I	IonInjectionTime';
len0 = length(keyword0);
len1 = length(keyword1);
len2 = length(keyword2);
len3 = length(keyword3);
len4 = length(keyword4);
%% get the MS1 info
% get the datatype
str=fgets(fid);
while feof(fid)==0 && 0==strcmp( str(1:len0),keyword0 )% If not at the end of the file, keep finding DataType
    str=fgets(fid);
end
MS1_datamode = str(len0+2);

if 1==strcmp('P',MS1_datamode)% DataType can be Centroid starting with C, or Profile starting with P
    fprintf(1,'MS1 is profile mode, convert to centroid mode first!\n');
    fclose(fid);
    return;
end

% for progress
ct_prt = 0;% length of output characters
fprintf(1,'MS1 scans: ');

% start to process
str = fgets(fid);
while 0==feof(fid)
    % If the symbol at the beginning of the spectrum (letter S) is found, start a round of reading, otherwise continue to read the next line, until the end of the file
    if 1==strcmp( str(1:len1),keyword1 )
        % progress
        fno = fno + 1;
        fprintf(repmat('\b',[1,ct_prt]));% Use backspace to go back to the existing output
        ct_prt = fprintf('%i',fno);% Output which one is being processed, progress prompt

        % 1.get the MS1 info
        % MS1 scan
        scan_no = str2num(str(len1+2:end));%#ok
        scan_no = scan_no(1);

        % RT
        str = fgets(fid);
        while feof(fid)==0 && 0==strcmp( str(1:len2),keyword2 )
            str=fgets(fid);
        end
        rt_no = str2double(str(len2+2:end));
        

        % IonInjectionTime
        str = fgets(fid);
        MS1_injecTime = str2double(str(len4+2:end));
        % InstrumentType
        str = fgets(fid);
		%
        if 1==fno
            MS1Type = str(len3+2:len3+5);
%             if 1==strcmp(MS1Type,'ITMS')
%                 fclose(fid);
%                 disp([ms1_fullfile,': low resolution MS1!']);
%                 return;
%             end;
        end
		%}

        % 2.read the MS1 data
        mz = zeros([1,maxpeaknum]);
        inten = zeros([1,maxpeaknum]);
        pnum = 0;
        while feof(fid)==0 && 0==strcmp( str(1:len1),keyword1 )% Keep reading until the current line starts with the letter S (indicating a new record)
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
        IX = find(inten>0);% Record the index of non-zero intensity
        mz = mz(IX);
        inten = inten(IX);

        % 3.judge whether to centroid MS1 scan

        % 4.save the MS1 info and peaks
        if 1==isempty(mz) || mz(end)-mz(1)<10
            fno = fno - 1;
            continue;
        end
        if length(inten)>100
            baseline = GetBaseline(inten);% After logarithmic transformation, create a histogram and take the center of the highest bar as the baseline
            IX = find(inten>=baseline);% Keep those with intensity higher than the baseline, which is a filtering process to remove low-abundance peaks
            mz = mz(IX);
            inten = inten(IX);
        else
            baseline = 0;
        end
        npk = length(mz);
        
        % scan_no is the spectrum starting positions in pe, rt_no is the retention time, npk is the number of peaks obtained after filtering low-abundance peaks,
        % baseline is the baseline used for filtering peaks, MS1_injecTimem is the IonInjectionTime in the MS1 file.
        MS1_index(fno,1:5) = [scan_no rt_no npk baseline MS1_injecTime];
        % pkno is the total number of peaks already recorded, the second parameter 1 is m/z, 2 is intensity
        MS1_peaks(pkno+1:pkno+npk,1) = mz;
        MS1_peaks(pkno+1:pkno+npk,2) = inten;
        pkno = pkno + npk;
    else
        str = fgets(fid);
    end
end
fclose(fid);% close the TXT file
fprintf(repmat('\b',[1,ct_prt]));
fprintf('%i',fno);
fprintf(1,'\n');

%% save the MS1 info
% filter the empty values
if fno<maxMS1num
    IX = 1:fno;
    MS1_index = MS1_index(IX,:);
end
tmp = MS1_index(1:fno,3);
MS1_index(1:fno,3) = cumsum(tmp) + 1;% Since we want to keep the row number of the next peak, why not keep the pkno variable from the beginning?
% if MS1_index(fno,2)>1000% 20*60=1200
%     MS1_index(1:fno,2) = MS1_index(1:fno,2)/60;% %If the retention time is too long, use minutes as the unit?
% end
MS1_index(1:fno,2) = MS1_index(1:fno,2)/60; % the retention time is saved in the unit of [minute]

if pkno<totalpeaknum
    IX = 1:pkno;
    MS1_peaks = MS1_peaks(IX,:);
end

% save the results
save(MS1_scanfile,'MS1_index','MS1Type');
save(MS1_peakfile,'MS1_peaks');

success = 1;
end



function baseline = GetBaseline(inten)
% Function: Get the baseline intensity for discarding peaks
% Input: intensity vector
% Output: baseline
% Method: Take logarithm, divide into steps, create histogram, take the center of the highest bar, restore baseline

loginten = log10(inten);
t = min(loginten):0.08:max(loginten);
[n,xout] = hist(loginten,t);

[tmp,idx] = max(n);%#ok
baseline = 10^xout(idx);
end