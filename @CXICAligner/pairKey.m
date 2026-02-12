function key = pairKey(~, ref_raw, target_raw)
% Build a canonical key for a ref/target run pair.
% Input:
%   ref_raw (1 x 1 char/string)
%       Reference raw name
%   target_raw (1 x 1 char/string)
%       Target raw name
% Output:
%   key (1 x 1 char)
%       Pair key string

key = [ref_raw, '->', target_raw];
end
