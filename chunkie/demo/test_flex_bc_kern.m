phi = 0.3;
src = []; src.r = [0;0]; src.n = [cos(phi);sin(phi)]; src.d = [cos(phi+pi/2);sin(phi+pi/2)];

targ = []; targ.r = [1;1]; targ.n = [1;1]/sqrt(2); targ.d = [1;-1]/sqrt(2);
kappa = 0.282842712474619;

% targ = []; targ.r = [1;1]; targ.n = [1;0]; targ.d = [1;0];
% targ = []; targ.r = [1;1]; targ.n = [0;1]; targ.d = [0;-1];

targ.d2 = [0.1;0.3];

zk = 0.9; nu = 0.3;
% nu = 2;
h = 1e-3;

ifree = 1;

%%
targh_n = [];
targh_n.r = targ.r + h*[1,0,-1].*targ.n;
targh_n.n = targ.n + 0*[1,0,-1];
targh_n.d = targ.d + 0*[1,0,-1];

if ifree
ikern = @(s,t) chnk.flex2d.kern(zk, s, t, 'free_plate_eval', nu); 
ikern_1 = @(s,t) chnk.flex2d.kern(zk, s, t, 'free_plate_bc1', nu); 
else
ikern = @(s,t) chnk.flex2d.kern(zk, s, t, 'clamped_plate_eval', nu); 
ikern_1 = @(s,t) chnk.flex2d.kern(zk, s, t, 'clamped_plate_bc1', nu); 
end

% ikern = @(s,t) chnk.flex2d.kern(zk, s, t, 's', nu); 
% ikern_1 = @(s,t) chnk.flex2d.kern(zk, s, t, 'sprime', nu);

u1 = ikern_1(src,targ);

uh= ikern(src,targh_n); u1_fd = (uh(1,:)-uh(3,:))/2/h;
[u1;u1_fd;abs(u1_fd-u1)]
assert(norm(u1_fd-u1)<1e-6)

% return
%%
targh_lap = [];
targh_lap.r = targ.r + h*[1,0,-1,0,0].*targ.n+h*[0,0,0,1,-1].*targ.d;
targh_lap.n = targ.n + 0*[1,0,-1,0,0];
targh_lap.d = targ.d + 0*[1,0,-1,0,0];
targh_lap.d2 = targ.d2 + 0*[1,0,-1,0,0];

if ifree
ikern_2 = @(s,t) chnk.flex2d.kern(zk, s, t, 'free_plate_bc2', nu); 
else
ikern_2 = @(s,t) chnk.flex2d.kern(zk, s, t, 'clamped_plate_bc2', nu); 
end
% ikern_2 = @(s,t) chnk.flex2d.kern(zk, s, t, 'supported_plate_bcs', nu); 
u2 = ikern_2(src,targ);%u2 = u2(2);

uh= ikern(src,targh_lap); 
u2_fd = (1-nu)*(uh(1,:)-2*uh(2,:)+uh(3,:))/h/h + ...
    nu*(uh(1,:)-4*uh(2,:)+uh(3,:)+uh(4,:)+uh(5,:))/h/h;
[u2;u2_fd;abs(u2_fd-u2)]
assert(norm(u2_fd-u2)<1e-6)


%%
if ~ifree
targh_lap = [];
nsten = [2,1,0,-1,-2]; dsten = [0,0,0,0,0];
targh_lap.r = targ.r + h*nsten.*targ.n+h*dsten.*targ.d;
targh_lap.n = targ.n + 0*nsten;
targh_lap.d = targ.d + 0*nsten;

ikern_3 = @(s,t) chnk.flex2d.kern(zk, s, t, 'clamped_plate_bc3', nu); 
% ikern_2 = @(s,t) chnk.flex2d.kern(zk, s, t, 'supported_plate_bcs', nu); 
u3 = ikern_3(src,targ);%u2 = u2(2);

uh= ikern(src,targh_lap); 
u3_fd_nnn = (0.5*uh(1,:)-uh(2,:)+uh(4,:)-0.5*uh(5,:))/h^3;

targh_lap.r = targ.r + h*nsten.*targ.n+h*targ.d;
uh_a= ikern(src,targh_lap);
u3_fd_nnt_a = (uh_a(2,:)-2*uh_a(3,:)+uh_a(4,:))/h^2;

targh_lap.r = targ.r + h*nsten.*targ.n-h*targ.d;
uh_b= ikern(src,targh_lap);
u3_fd_nnt_b = (uh_b(2,:)-2*uh_b(3,:)+uh_b(4,:))/h^2;

u3_fd_nnt = (u3_fd_nnt_a - u3_fd_nnt_b)/2/h;

u3_fd_nn = (uh(2,:)-2*uh(3,:)+uh(4,:))/h^2;

targh_lap.r = targ.r +h*nsten.*targ.d;
uh= ikern(src,targh_lap);
u3_fd_tt = (uh(2,:)-2*uh(3,:)+uh(4,:))/h^2;

u3_fd = u3_fd_nnn + (2-nu)*u3_fd_nnt + (1-nu)*kappa*(u3_fd_tt-u3_fd_nn);
% u3_fd = (u3_fd_nn);
[u3;u3_fd;abs(u3_fd-u3)]
assert(norm(u3_fd-u3)<1e-4)
% %%
% f = @(x) exp(2*x);
% fh = f(h*nsten);
% (0.5*fh(1)-fh(2)+fh(4)-0.5*fh(5))/h^3 - 2^3*f(0)
end