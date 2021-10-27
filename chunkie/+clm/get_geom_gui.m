function [chnkr,clmparams] = get_geom_gui(icase,opts)
clmparams = clm.setup(icase,opts);

if isfield(clmparams,'cpars')
  cpars = clmparams.cpars;
end
if isfield(clmparams,'cparams')
  cparams = clmparams.cparams;
end
if isfield(clmparams,'ncurve')
  ncurve = clmparams.ncurve;
end

if isfield(clmparams, 'nch')
  nch = clmparams.nch;
end

if isfield(clmparams,'k')
  k = clmparams.k;
end
if isfield(clmparams,'c')
  c = clmparams.c;
end
if isfield(clmparams,'coef')
  coef = clmparams.coef;
end

chnkr(1,ncurve) = chunker();

% define functions for curves
fcurve = cell(1,ncurve);
for icurve=1:ncurve
  fcurve{icurve} = @(t) clm.funcurve(t,icurve,cpars{icurve},icase);
end

% number of Gauss-Legendre nodes on each chunk
ngl = 16;
pref = []; 
pref.k = ngl;
disp(['Total number of unknowns = ',num2str(sum(nch)*ngl*2)])

% discretize the boundary

for icurve=1:ncurve
  chnkr(icurve) = chunkerfuncuni(fcurve{icurve},nch(icurve),cparams{icurve},pref);
  chnkr(icurve) = chnkr(icurve).makedatarows(4);
  chnkr(icurve).data(1,:,:) = k(c(1,icurve));
  chnkr(icurve).data(2,:,:) = k(c(2,icurve));
  chnkr(icurve).data(3,:,:) = coef(c(1,icurve));
  chnkr(icurve).data(4,:,:) = coef(c(1,icurve));
end

end