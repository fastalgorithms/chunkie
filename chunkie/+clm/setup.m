function [clmparams] = setup(icase,opts)

if icase == 2
  [clmparams] = clm.setup2();
elseif icase == 4
  [clmparams] = clm.setup4();
elseif icase == 6
  [clmparams] = clm.setup6(opts);
end
