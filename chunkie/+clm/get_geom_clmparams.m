function [chnkr] = get_geom_clmparams(clmparams)
ncurve = clmparams.ncurve;
chnkr(1,ncurve) = chunker();

% number of Gauss-Legendre nodes on each chunk
ngl = clmparams.ngl;
pref = []; 
pref.k = ngl;
disp(['Total number of unknowns = ',num2str(sum(clmparams.nch)*ngl*2)])

% discretize the boundary

for icurve=1:ncurve
  chnkr(icurve) = chunkerfuncuni(clmparams.fcurve{icurve},clmparams.nch(icurve),clmparams.cparams{icurve},pref);
  chnkr(icurve) = chnkr(icurve).makedatarows(4);
  chnkr(icurve).data(1,:,:) = clmparams.k(clmparams.c(1,icurve));
  chnkr(icurve).data(2,:,:) = clmparams.k(clmparams.c(2,icurve));
  chnkr(icurve).data(3,:,:) = clmparams.coef(clmparams.c(1,icurve));
  chnkr(icurve).data(4,:,:) = clmparams.coef(clmparams.c(1,icurve));
end

end