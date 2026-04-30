classdef CMS2ResultIO
    % Read/write helper for MS2 single-spectrum processing result

    methods (Static)
        keys = CMS2ResultFieldKeys()
        resultObj = read(msms_res_path)
        write(result, path, include_vif)
        cStr = formatImpStrings(massArrangement, fixedPosMod, dictVariMod, inxSites, pepSeq)
    end
end
