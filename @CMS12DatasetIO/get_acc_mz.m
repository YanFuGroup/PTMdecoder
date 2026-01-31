function acc_mz = get_acc_mz(obj,cen_mz,cur_mz,cur_chg)
% Calculate a more accurate mass-to-charge ratio
% Input:
%   obj (CMS12DatasetIO)
%       dataset IO instance
%   cen_mz (1 x 1 double) m/z
%       ActivationCenter
%   cur_mz (1 x 1 double) m/z
%       current m/z (original MS1 peak m/z)
%   cur_chg (1 x 1 double/int)
%       current charge
% Output:
%   acc_mz (1 x 1 double) m/z
%       peak closest to ActivationCenter; if none, use ActivationCenter

if obj.m_ms1_tolerance.isppm
    ptol = obj.m_ms1_tolerance.value*cen_mz*1e-6;
else
    ptol = obj.m_ms1_tolerance.value;
end
sets = [-1 0 1 2 3 4 5 6 7 8];
mzs = cur_mz + sets*CConstant.unitdiff/cur_chg;% m/z of various isotope peaks
ix = find(abs(mzs-cen_mz)<ptol*cen_mz*1e-6);
if 0==isempty(ix)
    [tmp,i] = min(abs( mzs(ix)-cen_mz ));%#ok
    acc_mz = mzs(ix(i));
else
    acc_mz = cen_mz;
end
end