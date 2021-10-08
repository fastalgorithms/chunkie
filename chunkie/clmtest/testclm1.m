function testclm1
%% This program solves the following layered medium problem
% legacy code, superseded by testclm2. Do not use!!!
%
%
%              Omega_1
%                 ___  ________
%                /   \/        \ 
%               /               \
%  -------------     Omega_3     ---------------------
%               \               /
%                \_____/\______/
%                  
%              Omega_2
%
%
%
%  \Delta u_i + k_i^2 u_i = 0, i=1,2,3
% 
%  u_i^{in}-u_i^{ex} = f_i
%
%  c_i^{in}du_i^{in}/dn - c_i^{ex}du_i^{ex}/dn = g_i
%
%  where "in" denotes interior, "ex" denotes exterior.
%
%  The unit normal vector is ALWAYS outward w.r.t. the interior domain.
%
%
%  The representations:
%
%  u_i = sum_{j\in B_i} u_{i,j},
%
%  u_{i,j} = c_j^{in/ex} D_{i,\Gamma_j}[\mu_j] + S_{i,\Gamma_j}[\rho_j],
%
%  where B_i contains all curves that form the boundary of the ith domain,
%  j is the curve index, c_j^{in/ex} = c_j^{in} if the ith domain is the 
%  exterior domain w.r.t. the jth curve, and c_j^{in/ex} = c_j^{ex} if the 
%  ith domain is the interior domain w.r.t. the jth curve. 
%
%
%
%  SKIEs:
%
%  (c_i^{in}+c_i^{ex})/2 \mu_i + (c_i^{in} D_i^{ex}[\mu_i]-c_i^{ex} D_i^{in}[\mu_i])
%     + (S_i^{ex}[\rho_i]-S_i^{in}[\rho_i])
%     - \sum_{j\ne i, j\in clist{c(1,i)}} u_{i,j}
%     + \sum_{j\ne i, j\in clist{c(2,i)}} u_{i,j}
%  = -f_i
%
%
%  (c_i^{in}+c_i^{ex})/2 \rho_i - c_i^{in} * c_i^{ex}(D'_i^{ex}[\mu_i]- D'_i^{in}[\mu_i])
%     - (c_i^{ex} S'_i^{ex}[\rho_i]- c_i^{in} S'_i^{in}[\rho_i])
%     + c_i^{in} \sum_{j\ne i, j\in clist{c(1,i)}} du_{i,j}/dn  
%     - c_i^{ex} \sum_{j\ne i, j\in clist{c(2,i)}} du_{i,j}/dn  
%  = g_i
%  
%
%
%
close all
format long e
format compact

addpaths_loc();

ndomain = 3; % number of domains
ncurve = 4; % number of curve segments
chnkr(1,ncurve) = chunker();

rn = zeros(ndomain,1); 
% rn(i) is the index of refraction of the ith domain
rn(1) = 1.0;
rn(2) = 3.1;
rn(3) = 2.2;

% k0 is the wave number in vacuum
k0 = 1;

% k(i) is the wave number for the ith domain
k = k0*rn;
k(1)=k(2);
% coefficients in the boundary conditions on normal derivatives
% (a) \Gamma_1 and \Gamma_2 -- Omega_1: interior; Omega_2: exterior
cin(1:2) = 1.1; cex(1:2) = 2.3;
% (b) \Gamma_3 -- Omega_1: exterior; Omega_3: interior
cin(3) = 1.4; cex(3) = 3.5;
% (c) \Gamma_4 -- Omega_2: exterior; Omega_3: interior
cin(4) = 2.5; cex(4) = 1.7;

cin = ones(1,ncurve);
cex = cin;

% domain indices for each curve
c = zeros(2,ncurve);
% interior domain for the ith curve
c(1,1:2) = 1; c(1,3:4) = 3;
% exterior domain for the ith curve 
c(2,1:2) = 2; c(2,3) = 1; c(2,4) = 2;

k1 = zeros(1,ncurve); % wave numbers for the interior domain
k2 = zeros(1,ncurve); % wave numbers for the exterior domain

for i=1:ncurve
  k1(i) = k(c(1,i));
  k2(i) = k(c(2,i));
end

% curve lists for each domain, positive means interior, 
% and negative means exterior.
clist = cell(1,ndomain);
clist{1} = [1,-3,2];
clist{2} = [-1,-4,-2];
clist{3} = [3,4];

% curve parameters
% center eye range on the real x-axis [a,b]
a = -1d1;
b = 1d1;
% length for complexification
C = 10;
% left flat curve [-L(1),a]
L(1) = C-a;
% right flat curve [b,L(2)]
L(2) = b+C;

% two circular arcs for the center eye for now
theta = zeros(1,ncurve);
% upper curve opening angle
theta(3) = pi;
% lower curve opening angle
theta(4) = pi;

% parameters for the complexification of left and right flat parts.
c1 = 8;
c2 = 3;

% number of Gauss-Legendre nodes on each chunk
ngl = 16;

pref = []; 
pref.k = ngl;

% number of chunks on each curve
nch = zeros(1,ncurve);

n0 = 8;

lambda = 2*pi/max(abs(k(1)),abs(k(2)));
nch(1) = round((a+L(1))/lambda) + n0;
nch(2) = round((L(2)-b)/lambda) + n0;

lambda = 2*pi/max(abs(k(1)),abs(k(3)));
nch(3) = round((b-a)/2/sin(theta(3)/2)/lambda) + n0;
lambda = 2*pi/max(abs(k(2)),abs(k(3)));
nch(4) = round((b-a)/2/sin(theta(4)/2)/lambda) + n0;

nch
% discretize domain

cparams = cell(1,ncurve);

for i=1:ncurve
    cparams{i}.ifclosed = false;
end

cparams{1}.ta = -L(1);
cparams{1}.tb = a;

cparams{2}.ta = b;
cparams{2}.tb = L(2);

for icurve=3:4
    cparams{icurve}.ta = -theta(icurve)/2;
    cparams{icurve}.tb = theta(icurve)/2;
end

cpars = zeros(3,ncurve);
cpars(:,1) = [L(1),c1,c2];
cpars(:,2) = [L(2),c1,c2];
cpars(:,3) = [a, b, theta(3)];
cpars(:,4) = [a, b, theta(4)];

start = tic; 
% complexified interface
for icurve=1:2
    chnkr(icurve) = chunkerfuncuni(@(t) clm.complexx1(t,icurve,cpars(:,icurve)),...
        nch(icurve),cparams{icurve},pref);
end

% center eye shape
for icurve=3:4
    chnkr(icurve) = chunkerfuncuni(@(t) clm.funcurve(t,icurve,cpars(:,icurve)), ...
        nch(icurve),cparams{icurve},pref);
end

t1 = toc(start);

fprintf('%5.2e s : time to build geo\n',t1)

for i=1:ncurve
    nch(i) = chnkr(i).nch;
end

% total number of discretization points
np = sum(nch)*pref.k

%[~,~,info] = sortinfo(chnkr);
%assert(info.ier == 0);

% plot geometry and data
figure(1)
clf
plot(chnkr,'r.')
axis equal

% figure(2)
% clf
% quiver(chnkr)
% axis equal
% pause

% now build the system matrix
M = zeros(2*np);

opts = [];
start = tic; 

% diagonal constant for each curve
alpha = zeros(1,ncurve);
for i=1:ncurve
  alpha(i) = 2/(cin(i)+cex(i));
end

for i=3:ncurve % source curve id
  %
  % define kernels
  sk1 =  @(s,t) chnk.helm2d.kern(k1(i),s,t,'s',1);
  sk2 =  @(s,t) chnk.helm2d.kern(k2(i),s,t,'s',1);

  dk1 =  @(s,t) chnk.helm2d.kern(k1(i),s,t,'d',1);
  dk2 =  @(s,t) chnk.helm2d.kern(k2(i),s,t,'d',1);

  spk1 = @(s,t) chnk.helm2d.kern(k1(i),s,t,'sprime',1);
  spk2 = @(s,t) chnk.helm2d.kern(k2(i),s,t,'sprime',1);

  dpk1 = @(s,t) chnk.helm2d.kern(k1(i),s,t,'dprime',1);
  dpk2 = @(s,t) chnk.helm2d.kern(k2(i),s,t,'dprime',1);
  
  indi1 = sum(nch(3:i-1))*ngl+(1:nch(i)*ngl);
  indi2 = indi1+np;
  
  for j=3:ncurve % target curve id
    if j==i
      % build matrices for 8 layer potentials
      S1  = chunkermat(chnkr(i),sk1,opts);
      S2  = chunkermat(chnkr(i),sk2,opts);

      D1  = chunkermat(chnkr(i),dk1,opts);
      D2  = chunkermat(chnkr(i),dk2,opts);

      S1p = chunkermat(chnkr(i),spk1,opts);
      S2p = chunkermat(chnkr(i),spk2,opts);

      D1p = chunkermat(chnkr(i),dpk1,opts);
      D2p = chunkermat(chnkr(i),dpk2,opts);
      
      M(indi1,indi1) =  alpha(i)*(cin(i)*D2-cex(i)*D1);
      M(indi1,indi2) =  alpha(i)*(S2-S1);
      M(indi2,indi1) = -alpha(i)*cex(i)*cin(i)*(D2p-D1p);
      M(indi2,indi2) = -alpha(i)*(cex(i)*S2p-cin(i)*S1p);
    else
      indj1 = sum(nch(3:j-1))*ngl+(1:nch(j)*ngl);
      indj2 = indj1+np;
%      if ismember(i,clist{c(1,j)}) || ismember(-i,clist{c(1,j)})
        if ismember(i,clist{c(1,j)})
          coeff = cex(i);
        else
          coeff = cin(i);
        end
        
        S1  = chunkermat_smooth(chnkr(i),chnkr(j),sk1);
        D1  = chunkermat_smooth(chnkr(i),chnkr(j),dk1);
        S1p = chunkermat_smooth(chnkr(i),chnkr(j),spk1);
        D1p = chunkermat_smooth(chnkr(i),chnkr(j),dpk1);       
 
        M(indj1,indi1) = M(indj1,indi1) - alpha(j)*coeff*D1;
        M(indj1,indi2) = M(indj1,indi2) - alpha(j)*S1;
        M(indj2,indi1) = M(indj2,indi1) + alpha(j)*cin(j)*coeff*D1p;
        M(indj2,indi2) = M(indj2,indi2) + alpha(j)*cin(j)*S1p;
 %     end
      
%      if ismember(i,clist{c(2,j)}) || ismember(-i,clist{c(2,j)})
        if ismember(i,clist{c(2,j)})
          coeff = cin(i);
        else
          coeff = cex(i);
        end
        
        S2  = chunkermat_smooth(chnkr(i),chnkr(j),sk2);
        D2  = chunkermat_smooth(chnkr(i),chnkr(j),dk2);
        S2p = chunkermat_smooth(chnkr(i),chnkr(j),spk2);
        D2p = chunkermat_smooth(chnkr(i),chnkr(j),dpk2);       
 
        M(indj1,indi1) = M(indj1,indi1) + alpha(j)*coeff*D2;
        M(indj1,indi2) = M(indj1,indi2) + alpha(j)*S2;
        M(indj2,indi1) = M(indj2,indi1) - alpha(j)*cex(j)*coeff*D2p;
        M(indj2,indi2) = M(indj2,indi2) - alpha(j)*cex(j)*S2p; 
%      end
    end
  end
end

% add the identity matrix
M = M + eye(2*np);

t1 = toc(start);

fprintf('%5.2e s : time to assemble matrix\n',t1)


% src(:,~i) are sources for the ith domain
src = [-1, 1.3, 0; max(chnkr(3).r(2,:),[],'all')+1,...
  min(chnkr(4).r(2,:),[],'all')-2.2, 0];
src
hold on;plot(src(1,:),src(2,:),'r*')

% construct artificial boundary data for testing purpose
rhs = zeros(2*np,1);

for i=3:ncurve
  din = c(1,i);
  dex = c(2,i);

  ind1 = sum(nch(3:i-1))*ngl+(1:nch(i)*ngl);
  ind2 = ind1+np;
  
  targnorm = chnk.normal2d(chnkr(i));
  nx = targnorm(1,:); nx=nx.';
  ny = targnorm(2,:); ny=ny.';
  
  for j=1:ndomain
    if j~=din
      [u1,gradu1]=chnk.helm2d.green(k1(i),src(:,j),chnkr(i).r(:,:));
      du1dn = gradu1(:,1).*nx + gradu1(:,2).*ny;

      rhs(ind1) = rhs(ind1) - alpha(i)*u1;
      rhs(ind2) = rhs(ind2) + alpha(i)*cin(i)*du1dn;
    %end
    
    %if j~=dex
    else
      [u2,gradu2]=chnk.helm2d.green(k2(i),src(:,j),chnkr(i).r(:,:));
      du2dn = gradu2(:,1).*nx + gradu2(:,2).*ny;

      rhs(ind1) = rhs(ind1) + alpha(i)*u2;
      rhs(ind2) = rhs(ind2) - alpha(i)*cex(i)*du2dn;
    end
  end
end

rhs = rhs(:);

% solve the linear system using gmres
start = tic; sol = gmres(M,rhs,[],1e-13,100); t1 = toc(start);

fprintf('%5.2e s : time for dense gmres\n',t1)

% compute the solution at one target point in each domain
targ1 = [-10; 1000]; % target point in domain #1
targ3 = [0; 0]; % target point in domain #2
targ2 = [2; -100]; % target point in domain #3

targ = [targ1, targ2, targ3]

sol1 = sol(1:np); % double layer density
sol2 = sol(np+1:end); % single layer density

uexact = zeros(ndomain,1);
ucomp = zeros(ndomain,1);

% compute the exact solution
for i=1:ndomain
  for j=1:ndomain
    %if j~=i
    if (i==3 && j~=3) || (i~=3 && j==3)
      uexact(i) = uexact(i) + chnk.helm2d.green(k(i),src(:,j),targ(:,i));
    end
  end
end

% compute the numerical solution
for i=1:ndomain
  skern =  @(s,t) chnk.helm2d.kern(k(i),s,t,'s',1);
  dkern =  @(s,t) chnk.helm2d.kern(k(i),s,t,'d',1);

  %for jcurve=1:length(clist{i})
   for j=3:ncurve 
    %j = clist{i}(jcurve);

    if j>0 % interior domain w.r.t. the jth curve
      coeff = cex(j);
    elseif j<0 % exterior domain w.r.t the jth curve
      j = -j;
      coeff = cin(j);
    end
   
    ind = sum(nch(3:j-1))*ngl+(1:nch(j)*ngl);
    
    dlp = chunkerkerneval(chnkr(j),dkern,sol1(ind),targ(:,i));
    slp = chunkerkerneval(chnkr(j),skern,sol2(ind),targ(:,i));
%     if i==1 && j==3
%       coeff=-coeff;
%     end
    
    ucomp(i) = ucomp(i) + coeff*dlp + slp;
  end
end

uerror = abs(ucomp-uexact)./abs(uexact);

[ucomp, uexact, uerror]

