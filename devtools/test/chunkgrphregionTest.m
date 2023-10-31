%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   .  .  .  builds a simple pentagonal chunkergraph 
%            and tests the interior Helmholtz Dirichlet problem
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%clear all

addpaths_loc();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   .  .  .  builds a simple pentagonal chunkergraph 
%            and tests the interior Helmholtz Dirichlet problem
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a = -1.0;
b = 1.0;

ncurve = 1;

% set angles for top and bottom curve
theta = [pi/2.2 pi/2.4];

% set parametrization start and end points for chunkie objects
tab = [-theta(1)/2 -theta(2)/2; theta(1)/2 theta(2)/2];

% set parameters for curve parametrization
verta = [-0.5,-0.5,0.5,0.5;-0.5,0.5,0.5,-0.5];
vertb = verta*2;
vertc = verta*3;
vertd = vertc;
vertc(1,:) = vertc(1,:)+12;

evert = [1.1;0];
verts = [vertc,vertd,vertb,verta,evert];
v2e = [1,-1,0,0;0,1,-1,0;0,0,1,-1;-1,0,0,1];
edge2verts = zeros(18,size(verts,2));
edge2verts(1:16,1:16) = kron(eye(4),v2e);
edge2verts(17,11) = 1;
edge2verts(17,17)=-1;
edge2verts(18,12) = 1;
edge2verts(18,17)= -1;

fchnks = [];


prefs      = [];
% prefs.chsmall = 1d-4;
[cgrph] = chunkgraph(verts,edge2verts,fchnks,prefs);

vstruc = procverts(cgrph);
cgrph = balance(cgrph);

% rgns = cgrph.regions;
% rgn1 = rgns{1};
% rgn2 = rgns{2};
% 
% r2ub = rgn2{1};
% 
% e2 = rgn2{1}{1}(1);
% v2 = find(cgrph.edge2verts(:,abs(e2))==1);
% v2 = verts(:,v2);
% 
% irgn = 0;
% for ii=1:numel(rgn1)
%     nin = pointinregion(cgrph,rgn1{ii},v2);
%     if (nin > 0 && mod(nin,2)==1)
%         if (irgn ~= 0)
%             display("Warning: an unsupported geometry error has occurred");
%         end
%         irgn = ii;
%     end
% end
% 
% if (irgn == 0)
%     
% end    
% 
% rgnout = rgns{1};
% if (numel(rgns)>1)
%     rgn2 = rgns{2};
%     [rgnout] = mergeregions(cgrph,rgnout,rgn2);
%     for ii=3:numel(rgns)
%         rgn2 = rgns{ii};
%         [rgnout] = mergeregions(cgrph,rgnout,rgn2);
%     end
% end    
%%%% [in] = inpolygon(xquery,yquery,xverts,yverts)



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
