chunkgrphrcipTransmissionTest0();


function chunkgrphrcipTransmissionTest0()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   .  .  .  builds a simple pentagonal chunkergraph 
%            and tests the Helmholtz transmission problem
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a = -1.0;
b = 1.0;

ncurve = 1;

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


verts = exp(1i*2*pi*(0:4)/5);
verts = [real(verts);imag(verts)];

[~, nv] = size(verts);

edgesendverts = [1:nv; [2:nv 1]];
[~, ncurve] = size(edgesendverts);
amp = 0.5;
frq = 6;

fchnks    = cell(1,ncurve);
cparams = cell(ncurve,1);
for icurve = 1:ncurve
    fchnks{icurve} = @(t) circulararc(t,cpars{1});
    fchnks{icurve} = @(t) sinearc(t,amp,frq);
    cparams{icurve}.ta = 0;
    cparams{icurve}.tb = 1;
end

[cgrph] = chunkgraph(verts, edgesendverts, fchnks, cparams);

plot(cgrph); hold on;
quiver(cgrph);

nregions = 2;
ks = [1.1;2.1]*10;
coefs = [1.0;1.0];
cs(1,1:ncurve) = 1;
cs(2,1:ncurve) = 2;
opts = [];
opts.bdry_data_type = 'point sources';

sources = cell(1,2);
sources{1} = [0.0;0.1];
sources{2} = [3.0;-3.2];
charges{1} = (1.2+1j)*10;
charges{2} = (1+0j)*10;

opts.sources = cell(1,2);
opts.sources{1} = sources{1};
opts.sources{2} = sources{2};
opts.charges = cell(1,2);
opts.charges{1} = charges{1};
opts.charges{2} = charges{2};

[kerns, bdry_data, kerns_eval] = chnk.helm2d.transmission_helper(cgrph, ...
                                   ks, cs, coefs, opts);


opts = [];
opts.nonsmoothonly = false;
opts.rcip = true;
[sysmat] = chunkermat(cgrph, kerns, opts);
sysmat = sysmat + eye(size(sysmat,2));


tic, dens = sysmat\bdry_data; toc;

%% Postprocessing

% generate some targets...

xs = -2:0.01:2;
ys = -2:0.01:2;
[X,Y] = meshgrid(xs,ys);
targs = [X(:).';Y(:).'];

opts_eval = [];
opts_eval.forcesmooth = true;
[utarg, targdomain] = chunkerkerneval(cgrph, kerns_eval, dens, targs, ...
      opts_eval);

[~, nt] = size(targs);
true_sol = zeros(nt,1);

iind1 = (targdomain == 1);
t.r = targs(:,iind1);
s = [];
s.r = sources{1};
true_sol(iind1) = charges{1}*chnk.helm2d.kern(ks(1), s, t, 's');


iind2 = (targdomain == 2);
t.r = targs(:,iind2);
s.r = sources{2};
true_sol(iind2) = charges{2}*chnk.helm2d.kern(ks(2),s,t,'s');


uerr = utarg - true_sol;
uerr = reshape(uerr,size(X));
figure
h = pcolor(X,Y,log10(abs(uerr)));
set(h,'EdgeColor','None'); hold on;
plot(cgrph,'w-','LineWidth',2);
caxis([-16,0])
colorbar

true_sol = reshape(true_sol,size(X));
figure
h = pcolor(X,Y,imag(true_sol));
set(h,'EdgeColor','None'); hold on;
plot(cgrph,'w-','LineWidth',2);
colorbar





end


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

function [r,d,d2] = sinearc(t,amp,frq)
xs = t;
ys = amp*sin(frq*t);
xp = ones(size(t));
yp = amp*frq*cos(frq*t);
xpp = zeros(size(t));
ypp = -frq*frq*amp*sin(t);

r = [(xs(:)).'; (ys(:)).'];
d = [(xp(:)).'; (yp(:)).'];
d2 = [(xpp(:)).'; (ypp(:)).'];
end
