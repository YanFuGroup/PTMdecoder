function obj = read(path)
% Read a report file into CIMPQuantReport
fin = fopen(path, 'r');
if fin < 0
    error(['Cannot open the report file "', path, '".']);
end

blocks = CIMPQuantBlock.empty(0,1);
current_protein_name_pos = {};
current_records = CIMPQuantRecord.empty(0,1);
current_record = CIMPQuantRecord.empty(0,1);
has_record = false;
has_block = false;

% Skip header lines
fgetl(fin);
fgetl(fin);
fgetl(fin);

while ~feof(fin)
    strline = fgetl(fin);
    if ~ischar(strline) || isempty(strline)
        continue;
    end

    if strline(1) == '@'
        if ~has_record
            continue;
        end
        peak = CIMPQuantResultIO.parse_rt_peak_line(strline);
        current_record.rt_peaks(end+1,1) = peak; %#ok<AGROW>
    elseif strline(1) == '*'
        if has_record
            current_records(end+1,1) = current_record; %#ok<AGROW>
        end
        current_record = CIMPQuantResultIO.parse_record_line(strline);
        has_record = true;
    else
        if has_record
            current_records(end+1,1) = current_record; %#ok<AGROW>
            has_record = false;
        end
        if has_block
            blocks(end+1,1) = CIMPQuantBlock(current_protein_name_pos, current_records); %#ok<AGROW>
        end
        current_protein_name_pos = CIMPQuantResultIO.parse_protein_line(strline);
        current_records = CIMPQuantRecord.empty(0,1);
        has_block = true;
    end
end

if has_record
    current_records(end+1,1) = current_record; %#ok<AGROW>
end
if has_block
    blocks(end+1,1) = CIMPQuantBlock(current_protein_name_pos, current_records); %#ok<AGROW>
end
fclose(fin);

obj = CIMPQuantReport(blocks);
end
