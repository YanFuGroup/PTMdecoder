classdef CMS2ResultIO
    % Read/write helper for MS2 single-spectrum processing result

    methods (Static)
        resultObj = read(msms_res_path)
        write(result, path)
        cStr = formatImpStrings(massArrangement, fixedPosMod, dictVariMod, inxSites, pepSeq)
    end
end
