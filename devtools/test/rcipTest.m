
%RCIPTEST
%
% This file tests the rcip routines for solving the exterior dirichlet 
% problem on a domain defined by two arcs of circles meeting at two vertices

clearvars; close all;
addpaths_loc();

ncurve = 2;
chnkr(1,ncurve) = chunker();

% set wave number
zk = 1.1;




nch = 5*ones(1,ncurve);

a = -1.0;
b = 1.0;

% set angles for top and bottom curve
theta = [pi/2.2 pi/2.4];

% set parametrization start and end points for chunkie objects
tab = [-theta(1)/2 -theta(2)/2; theta(1)/2 theta(2)/2];

% set parameters for curve parametrization
vert = [a b; 0 0];
cpars = cell(1,ncurve);
cpars{1}.v0 = vert(:,1); cpars{1}.v1 = vert(:,2);
cpars{1}.theta = theta(1); cpars{1}.ifconvex = 0; cpars{1}.islocal = -1;

cpars{2}.v0 = vert(:,1); cpars{2}.v1 = vert(:,2);
cpars{2}.theta = theta(2); cpars{2}.ifconvex = 2; cpars{2}.islocal = -1;

% number of gauss legendre nodes
ngl = 16;
pref = [];
pref.k = ngl;

cparams = cell(1,ncurve);
for i=1:ncurve
    cparams{i}.ta = tab(1,i); 
    cparams{i}.tb = tab(2,i);
    cparams{i}.ifclosed = false;
end

inds = [0, cumsum(nch)];

[glxs,glwts,glu,glv] = lege.exps(ngl);

fcurve = cell(1,ncurve);
figure(1)
clf; hold on;
for icurve = 1:ncurve
    fcurve{icurve} = @(t) circulararc(t,cpars{icurve});
    chnkr(icurve) = chunkerfuncuni(fcurve{icurve},nch(icurve),cparams{icurve},pref);
    plot(chnkr(icurve)); quiver(chnkr(icurve));
end

iedgechunks = [1 2; 1 chnkr(2).nch];
isstart = [1 0];
nedge = 2;
ndim=1;
[Pbc,PWbc,starL,circL,starS,circS,ilist] = chnk.rcip.setup(ngl,ndim, ...
      nedge,isstart);

fkern = @(s,t) -2*chnk.helm2d.kern(zk,s,t,'D');

%%

% sources

ns = 10;
ts = 0.0+2*pi*rand(1,ns);
sources = [cos(ts); sin(ts)];
sources = 3.0*sources;
strengths = randn(ns,1);

% targets

nt = 20;
ts = 0.0+2*pi*rand(1,nt);
targets = [cos(ts);sin(ts)];
targets = 0.2*targets;

scatter(sources(1,:),sources(2,:),'o');
scatter(targets(1,:),targets(2,:),'x');
axis equal 




chnkrtotal = merge(chnkr);
fkern = @(s,t) -2*chnk.helm2d.kern(zk,s,t,'D');
np = chnkrtotal.k*chnkrtotal.nch;
start = tic; D = chunkermat(chnkrtotal,fkern);
t1 = toc(start);



kerns = @(s,t) chnk.helm2d.kern(zk,s,t,'s');

% eval u on bdry

targs = chnkrtotal.r; targs = reshape(targs,2,chnkrtotal.k*chnkrtotal.nch);
targstau = tangents(chnkrtotal); 
targstau = reshape(targstau,2,chnkrtotal.k*chnkrtotal.nch);

srcinfo = []; srcinfo.r = sources; 
targinfo = []; targinfo.r = targs;
kernmats = kerns(srcinfo,targinfo);
ubdry = kernmats*strengths;


% eval u at targets

targinfo = []; targinfo.r = targets;
kernmatstarg = kerns(srcinfo,targinfo);
utarg = kernmatstarg*strengths;



fprintf('%5.2e s : time to assemble matrix\n',t1)


sysmat = D;

RG = speye(np);
ncorner = 2;
corners= cell(1,ncorner);
R = cell(1,ncorner);


corners{1}.clist = [1,2];
corners{1}.isstart = [1,0];
corners{1}.nedge = 2;
corners{1}.iedgechunks = [1 2; 1 chnkr(2).nch];

corners{2}.clist = [1,2];
corners{2}.isstart = [0 1];
corners{2}.nedge = 2;
corners{2}.iedgechunks = [1 2; chnkr(1).nch 1];

ndim = 1;

nsub = 100;

opts = [];

for icorner=1:ncorner
    clist = corners{icorner}.clist;
    isstart = corners{icorner}.isstart;
    nedge = corners{icorner}.nedge;
    
    
    starind = zeros(1,2*ngl*ndim*nedge);
    for i=1:nedge
        i1 = (i-1)*2*ngl*ndim+1;
        i2 = i*2*ngl*ndim;
        if(isstart(i))
            
            starind(i1:i2) = inds(clist(i))*ngl*ndim+(1:2*ngl*ndim);
        else
            
            starind(i1:i2) = inds(clist(i)+1)*ngl*ndim-fliplr(0:2*ngl*ndim-1);
        end
    end
    
    
    [Pbc,PWbc,starL,circL,starS,circS,ilist] = chnk.rcip.setup(ngl,ndim, ...
      nedge,isstart);
    opts_use = [];
    
    iedgechunks = corners{icorner}.iedgechunks;
    tic; R{icorner} = chnk.rcip.Rcompchunk(chnkr,iedgechunks,fkern,ndim, ...
        Pbc,PWbc,nsub,starL,circL,starS,circS,ilist,... 
        glxs);
    toc
    
    sysmat(starind,starind) = inv(R{icorner}) - eye(2*ngl*nedge*ndim);
end
sysmat = sysmat + eye(np);

sol = gmres(sysmat,ubdry,np,eps*20,np);

opts.usesmooth=true;
opts.verb=false;
opts.quadkgparams = {'RelTol',1e-16,'AbsTol',1.0e-16};
start=tic; Dsol = chunkerkerneval(chnkrtotal,fkern,sol,targets,opts); 
t1 = toc(start);
fprintf('%5.2e s : time to eval at targs (slow, adaptive routine)\n',t1)

%

wchnkr = weights(chnkrtotal);

relerr = norm(utarg-Dsol,'fro')/(sqrt(chnkrtotal.nch)*norm(utarg,'fro'));
relerr2 = norm(utarg-Dsol,'inf')/dot(abs(sol(:)),wchnkr(:));
fprintf('relative frobenius error %5.2e\n',relerr);
fprintf('relative l_inf/l_1 error %5.2e\n',relerr2);


%%


%%
%----------------------------------------
%
%
% Auxiliary routines for generating boundary
%

function [r,d,d2] = circulararc(t,cpars)
%%circulararc
% return position, first and second derivatives of a circular arc
% that passes through points (x,y)=(a,0) and (x,y)=(b,0) with opening
% angle theta0.
%
% Inputs:
% t - paramter values (-theta0/2,theta0/2) to evaluate these quantities
%
% Outputs:
% r - coordinates
% d - first derivatives w.r.t. t
% d2 - second derivatives w.r.t. t

v0 = cpars.v0;
v1 = cpars.v1;
theta0 = cpars.theta;
ifconvex = cpars.ifconvex;

a = v0(1);
b = v1(1);

islocal = -1;
if isfield(cpars,'islocal')
  islocal = cpars.islocal; 
end

if islocal == -1
  cx = (a+b)/2;
  r0 = (b-a)/2/sin(theta0/2);

  if ifconvex 
    cy = v0(2)+(b-a)/2/tan(theta0/2);
    theta = 3*pi/2+t;
  else
    cy = v0(2)-(b-a)/2/tan(theta0/2);
    theta = pi/2+t;
  end


  xs = r0*cos(theta);
  ys = r0*sin(theta);

  xp = -ys;
  yp = xs;

  xpp = -xs;
  ypp = -ys;

  xs = cx + xs;
  ys = cy + ys;
  
else

  r0 = (b-a)/2/sin(theta0/2);  
  if ~ifconvex
    if islocal == 0
      cx = r0*cos(pi/2+theta0/2);
      sx = r0*sin(pi/2+theta0/2);
    elseif islocal == 1
      cx = r0*cos(pi/2-theta0/2);
      sx = r0*sin(pi/2-theta0/2);     
    end
  elseif ifconvex
    if islocal == 0
      cx = r0*cos(3*pi/2+theta0/2);
      sx = r0*sin(3*pi/2+theta0/2);
    elseif islocal == 1
      cx = r0*cos(3*pi/2-theta0/2);
      sx = r0*sin(3*pi/2-theta0/2);     
    end   
  end

%   t2 = t.*t;
%   t3 = t2.*t;
%   t4 = t3.*t;
%   t5 = t4.*t;
%   t6 = t5.*t;
%   t7 = t6.*t;
%   t8 = t7.*t;
%   t9 = t8.*t;
%   t10 = t9.*t;
%   
% 
%   ct = -t2/2 + t4/24 - t6/720 + t8/40320 - t10/3628800;
%   st =  t - t3/6 + t5/120 - t7/5040 + t9/362880;

% better to use sin(t) directly instead of its series expansion because
% 1. the evaluation of sin(t) is accurate at t=0;
% 2. one needs to determine how many terms are needed in series expansion,
% which looks simple but it could be either inefficient or inaccurate
% depending on the value of t; 3. if we write "if" statements, it's
% difficult to vectorize the evaluation.

% same reason that we should use -2*sin(t/2).^2 instead of its Taylor
% expansion, but of course one should NOT use cos(t)-1.
  ct = -2*sin(t/2).^2;
  st = sin(t);
  
  ctp = -st;
  stp = cos(t);
  
  ctpp = -stp;
  stpp = -st;
  
  xs = -sx*st + cx*ct;
  ys =  cx*st + sx*ct;
  
  xp = -sx*stp + cx*ctp;
  yp =  cx*stp + sx*ctp;

  xpp = -sx*stpp + cx*ctpp;
  ypp =  cx*stpp + sx*ctpp;
end    
    

r = [(xs(:)).'; (ys(:)).'];
d = [(xp(:)).'; (yp(:)).'];
d2 = [(xpp(:)).'; (ypp(:)).'];

end



