function resultObj = read(msms_res_path)
% Read from a msms result file
% Input:
%   msms_res_path (1 x 1 char/string)
%       path of the msms result file
% Output:
%   resultObj (CMS2Result)
%       CMS2Result object containing the parsed data

fieldKeys = CMS2ResultIO.CMS2ResultFieldKeys();

% Initialize result object
resultObj = CMS2Result();

fin = fopen(msms_res_path, 'r');
if fin < 0
    CLogger.error('[CMS2ResultIO:OpenFileFailed] Cannot open the msms level result:"%s"!', ...
        msms_res_path);
end
cleanup = onCleanup(@() fclose(fin));

file_total_length = dir(msms_res_path).bytes;
if file_total_length == 0
    CLogger.error('[CMS2ResultIO:EmptyFile] Warning: The file "%s" is empty!', ...
        msms_res_path);
end

% Read the file
while(~feof(fin))
    strline = fgetl(fin);
    if isempty(strline)
        continue;
    end
    tokens = strsplit(strline, '\t');
    tag = tokens{1};
    if strcmp(tag, 'P')
        % Record one peptide line (allow repeated P lines for the same sequence)
        if numel(tokens) < 2
            CLogger.error('[CMS2ResultIO:InvalidPeptideLine] Invalid peptide line: %s', strline);
        end
        resultObj.addOrSelectPeptide(tokens{2});
    elseif strcmp(tag, 'S')
        % Record one spectrum line
        if numel(tokens) < 3
            CLogger.error('[CMS2ResultIO:InvalidSpectrumLine] Invalid spectrum line: %s', strline);
        end

        spectrumMeta = struct('jaccard_stability', parseNamedNumericFieldFromTokensFromIdx( ...
            tokens, 4, fieldKeys.spectrum.jaccard, strline, ...
            'CMS2ResultIO:InvalidSpectrumNamedField', NaN), ...
            'vif_all_imp_max', parseNamedNumericFieldFromTokensFromIdx( ...
            tokens, 4, fieldKeys.spectrum.vifAll, strline, ...
            'CMS2ResultIO:InvalidSpectrumNamedField', NaN), ...
            'vif_reported_imp_max', parseNamedNumericFieldFromTokensFromIdx( ...
            tokens, 4, fieldKeys.spectrum.vifReported, strline, ...
            'CMS2ResultIO:InvalidSpectrumNamedField', NaN));

        resultObj.addSpectrum(tokens{2}, tokens{3}, spectrumMeta);
    else
        % Record one peptidoform line
        if numel(tokens) < 2
            CLogger.error('[CMS2ResultIO:InvalidPeptidoformLine] Invalid peptidoform line: %s', strline);
        end

        abundance = str2double(tokens{2});
        if isnan(abundance) && ~strcmpi(strtrim(tokens{2}), 'nan')
            CLogger.error('[CMS2ResultIO:InvalidPeptidoformLine] Invalid peptidoform line: %s', strline);
        end

        peptidoformMeta = struct('support_frequency', parseNamedNumericFieldFromTokensFromIdx( ...
            tokens, 3, fieldKeys.peptidoform.support, strline, ...
            'CMS2ResultIO:InvalidPeptidoformNamedField', NaN), ...
            'vif', parseNamedNumericFieldFromTokensFromIdx( ...
            tokens, 3, fieldKeys.peptidoform.vif, strline, ...
            'CMS2ResultIO:InvalidPeptidoformNamedField', NaN), ...
            'abundance_mad', parseNamedNumericFieldFromTokensFromIdx( ...
            tokens, 3, fieldKeys.peptidoform.mad, strline, ...
            'CMS2ResultIO:InvalidPeptidoformNamedField', NaN));

        resultObj.addPeptidoform(tokens{1}, abundance, peptidoformMeta);
    end
end

% Organize the data structure: compress buffers
resultObj.compress();
end


function value = parseNamedNumericFieldFromTokensFromIdx(tokens, startIdx, targetKey, rawLine, errorId, defaultValue)
% parseNamedNumericFieldFromTokensFromIdx - Parse a known named numeric field from tab-separated tokens
% Input:
%   tokens (1 x N cell)
%       tab-separated tokens from one raw line
%   startIdx (1 x 1 double)
%       first token index to scan for named fields
%   targetKey (1 x N char/string)
%       expected named-field key (for example: "jaccard", "support", "vif", or "mad")
%   rawLine (1 x N char/string)
%       original raw line, used in error reporting
%   errorId (1 x N char/string)
%       logical error id tag embedded in CLogger.error message
%   defaultValue (1 x 1 double)
%       default value returned when targetKey is missing
% Output:
%   value (1 x 1 double)
%       parsed numeric value for targetKey; defaultValue when not found
% Notes:
%   - Unknown named fields are ignored silently.
%   - Invalid numeric value for targetKey triggers CLogger.error.
value = defaultValue;
for t = startIdx:numel(tokens)
    token = strtrim(tokens{t});
    if isempty(token)
        continue;
    end

    kv = strsplit(token, '=');
    if numel(kv) ~= 2
        continue;
    end

    key = strtrim(kv{1});
    if strcmp(key, targetKey)
        value = str2double(strtrim(kv{2}));
        if isnan(value) && ~strcmpi(strtrim(kv{2}), 'nan')
            CLogger.error('[%s] Invalid named field value in line: %s', errorId, rawLine);
        end
    end
end
end