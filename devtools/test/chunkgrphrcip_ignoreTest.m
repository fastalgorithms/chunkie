chunkgrphrcip_ignoreTest0();


function chunkgrphrcip_ignoreTest0()
% chunkgrphrcip_ignoreTest verify that solution is accurate when
% artificial vertices are ignored



% set wave number
zk = 1.1;
ncurve = 3;

a = -1.0;
b = 1.0;

% set angles for top and bottom curve
theta = [pi/2.2 pi/2.4];


% set parametrization start and end points for chunkie objects
tab = [-theta(1)/2 -theta(2)/2 0; theta(1)/2 0 theta(2)/2];

% set parameters for curve parametrization
verts = [a 0 b; 0 (1/sin(theta(2)/2)-cot(theta(2)/2)) 0];
cpars = cell(1,ncurve);
cpars{1}.v0 = verts(:,1); cpars{1}.v1 = verts(:,3);
cpars{1}.theta = theta(1); cpars{1}.ifconvex = 0; cpars{1}.islocal = -1;

cpars{2}.v0 = verts(:,1); cpars{2}.v1 = verts(:,2);
cpars{2}.theta = theta(2)/2; cpars{2}.ifconvex = 2; cpars{2}.islocal = -1;

cpars{3}.v0 = verts(:,2); cpars{3}.v1 = verts(:,3);
cpars{3}.theta = theta(2)/2; cpars{3}.ifconvex = 2; cpars{3}.islocal = -1;
edgesendverts = [[1;3],[3;2],[2;1]];

% number of gauss legendre nodes
pref = [];
pref.k = 16;


cparams = cell(1,ncurve);
fchnks = cell(1,ncurve);
for i=1:ncurve
    cparams{i}.ta = tab(1,i); 
    cparams{i}.tb = tab(2,i);
    cparams{i}.ifclosed = false;

    fchnks{i} = @(t) circulararc(t,cpars{i});
end

cgrph = chunkgraph(verts,edgesendverts,fchnks,cparams,pref);

% sources

ns = 10;
ts = 0.0+2*pi*rand(1,ns);
sources = [cos(ts); sin(ts)];
sources = 3.0*sources;
strengths = randn(ns,1);

% targets

nt = 20;
ts = 0.0+2*pi*rand(1,nt);
targets = [2*cos(ts);sin(ts)];
targets = 0.2*targets;
% targets(:,1) = [0.95;0];
% targets(:,2) = [0,0.36];
% targets(:,3) = [-0.95;0];
% targets(:,2) = [0,0.36];

figure(1);clf
plot(cgrph)
daspect([1,1,1])
hold on
scatter(sources(1,:),sources(2,:),'o');
scatter(targets(1,:),targets(2,:),'x');
hold off

dkern = -2*kernel('helmholtz','d',zk);

% test with correct settings
opts = [];
opts.rcip = true;
opts.rcip_ignore = 2;
sys = chunkermat(cgrph,dkern,opts);
sys = eye(cgrph.npt) + sys;

skern = kernel('helmholtz','s',zk);
rhs = skern.eval(struct("r",sources),cgrph)*strengths;

sol = sys\rhs;

utrue = skern.eval(struct("r",sources),struct("r",targets))*strengths;
u = chunkerkerneval(cgrph,dkern,sol,targets);
assert(norm(u-utrue)<1e-10)

% test with incorrect settings
opts = [];
opts.rcip = true;
opts.rcip_ignore = 1;
sys = chunkermat(cgrph,dkern,opts);
sys = eye(cgrph.npt) + sys;

skern = kernel('helmholtz','s',zk);
rhs = skern.eval(struct("r",sources),cgrph)*strengths;

sol = sys\rhs;

utrue = skern.eval(struct("r",sources),struct("r",targets))*strengths;
u = chunkerkerneval(cgrph,dkern,sol,targets);
assert(norm(u-utrue)>1e-10)



% rmin = min(cgrph); rmax = max(cgrph);
% xl = rmax(1)-rmin(1);
% yl = rmax(2)-rmin(2);
% nplot = 400;
% xtarg = linspace(rmin(1),rmax(1),nplot); 
% ytarg = linspace(rmin(2),rmax(2),nplot);
% [xxtarg,yytarg] = meshgrid(xtarg,ytarg);
% targets = zeros(2,length(xxtarg(:)));
% targets(1,:) = xxtarg(:); targets(2,:) = yytarg(:);
% 
% %
% 
% start = tic; in = chunkerinterior(cgrph,{xtarg,ytarg}); t1 = toc(start);
% out = ~in;
% 
% fprintf('%5.2e s : time to find points in domain\n',t1)
% 
% % compute layer potential based on oversample boundary
% 
% start = tic;
% uscat = chunkerkerneval(cgrph,dkern,sol,targets(:,in)); t1 = toc(start);
% fprintf('%5.2e s : time for kernel eval (for plotting)\n',t1)
% 
% uin = skern.eval(struct("r",sources),struct("r",targets(:,in)))*strengths;
% utot = uscat(:)-uin(:);
% 
% 
% maxin = max(abs(uin(:)));
% maxsc = max(abs(uin(:)));
% maxtot = max(abs(uin(:)));
% 
% maxu = max(max(maxin,maxsc),maxtot);
% 
% figure(2)
% clf
% 
% zztarg = nan(size(xxtarg));
% zztarg(in) = utot;
% h=pcolor(xxtarg,yytarg,log10(abs(zztarg)));
% set(h,'EdgeColor','none')
% colorbar
% hold on
% plot(cgrph,'.k','LineWidth',2)
% axis equal tight
% set(gca, "box","off","Xtick",[],"Ytick",[]);
% 
% title('$\log_{10} error$','Interpreter','latex','FontSize',12)





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



