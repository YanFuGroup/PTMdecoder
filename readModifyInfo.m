% Create a dictionary mapping modification names to mono masses using the pFind studio style modification library file
% Input: strPathname (1 x N char/string) - the file name
% Output: mapModification (containers.Map) - a dictionary mapping modification names (with amino acid specificity) to mono mass offsets
function mapModification=readModifyInfo(strPathname)
% Compatibility wrapper; canonical implementation lives in CModificationRegistry.
mapModification = CModificationRegistry.loadMap(strPathname);
end

