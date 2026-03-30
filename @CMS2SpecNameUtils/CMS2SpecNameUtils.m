classdef CMS2SpecNameUtils
    % Utility helpers for MS2 spectrum-name parsing.

    methods (Static)
        function MS2ScanI = parseMS2ScanNumber(spectrum_name)
            % Parse scan number from spectrum name.
            % Inputs:
            %    spectrum_name (1 x 1 char/string)
            %        Spectrum name, either dotted form or raw scan number.
            % Outputs:
            %    MS2ScanI (1 x 1 double)
            %        Parsed MS2 scan number.

            spec_name = regexp(char(spectrum_name), '\.', 'split');
            if numel(spec_name) >= 2
                scan_str = spec_name{2};
            else
                scan_str = spec_name{1};
            end

            MS2ScanI = str2double(scan_str);
            if isnan(MS2ScanI)
                error('CMS2SpecNameUtils:InvalidSpectrumName', ...
                    'Cannot parse scan number from spectrum_name="%s".', char(spectrum_name));
            end
        end
    end
end
