function [chnkr,clmparams] = get_geom_gui(icase)
clmparams = clm.setup(icase);

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
end

end