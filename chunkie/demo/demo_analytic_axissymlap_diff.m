clearvars; clc;

iftorus = 0;
cparams = [];
cparams.eps     = 1.0e-10;
cparams.nover   = 1;

if ~iftorus  % sphere / hemisphere
    cparams.ta       = -pi/2;
    cparams.tb       =  pi/2;
    center           = [0; 2];
    cparams.ifclosed = false;
else          % torus
    cparams.ta       = 0;
    cparams.tb       = 2*pi;
    center           = [3; 0];
    cparams.ifclosed = true;
end

cparams.maxchunklen = 0.5;
radius = 1;
fcurve = @(t) radius*[cos(t(:).'); sin(t(:).')];

chnkrhalf = chunkerfunc(fcurve, cparams);
chnkrhalf = chnkrhalf.move(-center);

ndim  = 5;

ftrue  = @(s)  s.r(1,:).^2 - (ndim-1)*s.r(2,:).^2;

% First derivatives
dudr   = @(s)  2*s.r(1,:);
dudz   = @(s) -2*(ndim-1)*s.r(2,:);

% Second derivatives 
d2udrr  =  2;
d2udzz  = -2*(ndim-1);
d2udrz  =  0;

% Normal-derivative
dudn   = @(s)  dudr(s).*s.n(1,:) + dudz(s).*s.n(2,:);   % grad u . n

% Second normal derivative:  n^T H(u) n  
% H(u) = diag([2, -2*(ndim-1)])  is the Hessian
d2udn2 = @(s)  d2udrr*s.n(1,:).^2 + 2*d2udrz*s.n(1,:).*s.n(2,:) ...
             + d2udzz*s.n(2,:).^2;

% Mixed n_src / n_targ second derivative:  n_t^T H(u) n_s
d2udn_tdn_s = @(st, ss) ...
    d2udrr*st.n(1,:).'.*ss.n(1,:) + ...
    d2udrz*(st.n(2,:).'.*ss.n(1,:) + st.n(1,:).'.*ss.n(2,:)) + ...
    d2udzz*st.n(2,:).'.*ss.n(2,:);

% Solve interior Dirichlet problem
rhs    = 2*ftrue(chnkrhalf).';
kernd  = kernel('axissymlaplace', 'd',  ndim);
kerndp = kernel('axissymlaplace', 'dp', ndim);

mat  = chunkermat(chnkrhalf, 2*kernd);
A    = mat + eye(chnkrhalf.npt);
sig  = gmres(A, rhs, [], 1e-14, 200);
wts  = sig .* chnkrhalf.wts(:);   % weighted density
fprintf('Density solve residual: %e\n', norm(A*sig - rhs)/norm(rhs));

x1 = linspace(center(1),         center(1)+radius, 200);
x2 = linspace(center(2)-radius,  center(2)+radius, 200);
[xx, yy] = meshgrid(x1, x2);
in = chunkerinterior(chnkrhalf, {x1, x2});

t     = [];
t.r   = [xx(:).'; yy(:).'];
t.r   = t.r(:, in(:));
t.n   = repmat([0; 1], 1, size(t.r, 2));

%% check D
ucomp  = kernd.eval(chnkrhalf, t) * wts;
utrue  = ftrue(t).';
err_D  = log10(abs((ucomp - utrue) ./ max(abs(utrue))));
% fprintf('[D]   max log10 rel-err = %.2f\n', max(err_D));
% plot error
plotdata1 = nan(size(xx));
plotdata1(in) = err_D;
plot(chnkrhalf); hold on; plot(chnkrhalf,'bo'); quiver(chnkrhalf,'r');
h = pcolor(xx,yy,reshape(plotdata1,size(xx))); set(h,'EdgeColor','none');
title('err in D'); 
clim([-10,-1]);
colorbar; axis square;

%% check S and Sp
% Solve the interior Neumann problem u = S*sigma
%   (I + 2*S')*mu = 2*g   where g = du/dn on bdry
fprintf('\n--- Testing S and Sp ---\n');

% Neumann data
kersp  = kernel('axissymlaplace', 'sp', ndim);
kers = kernel('axissymlaplace', 's', ndim);
g_neu  = -2*dudn(chnkrhalf).';   % Neumann data on boundary
mat_sp  = chunkermat(chnkrhalf, -2*kersp);
A_sp    = eye(chnkrhalf.npt) + mat_sp;
mu      = gmres(A_sp, g_neu, [], 1e-13, 300);
wts_mu  = mu .* chnkrhalf.wts(:);
fprintf('Neumann solve residual: %e\n', norm(A_sp*mu - g_neu)/norm(g_neu));

% Verify single-layer recovers u
uS   = kers.eval(chnkrhalf,t) * wts_mu;
err_S = log10(abs((uS-utrue)./max(abs(utrue))));
%fprintf('[S]   max log10 rel-err = %.2f\n', max());

% plot error
plotdata1 = nan(size(xx));
plotdata1(in) = err_S;
plot(chnkrhalf); hold on; plot(chnkrhalf,'bo'); quiver(chnkrhalf,'r');
h = pcolor(xx,yy,reshape(plotdata1,size(xx))); set(h,'EdgeColor','none');
title('err in S'); 
clim([-10,-1]);
colorbar; axis square;


%% Test for Spp
kerspp = kernel('axissymlaplace', 'spp', ndim);
% Analytic solution  d_{n_t}^2 u(x)  with  n_t = [0;1]  => d^2u/dz^2 = -2*(ndim-1)
d2u_nt2_true = d2udn2(t).';   % constant = -2*(ndim-1) for n_t=[0,1]

uSpp  = kerspp.eval(chnkrhalf, t) * wts_mu;
err_Spp = log10(abs((uSpp - d2u_nt2_true) ./ max(abs(d2u_nt2_true))));
fprintf('[Spp]  max log10 rel-err = %.2f\n', max(err_Spp));

%% check Q
% u = S[mu]:
% Q[mu](x) = d_{n_y}^2 S[mu](x)
% adjoint of S'' in the L^2 sense.  
fprintf('\n--- Testing Q ---\n');

kernq = kernel('axissymlaplace', 'q', ndim);
uQ = kernq.eval(chnkrhalf, t) * wts_mu; 

% adjoint check:
% <Q[mu], phi>  =  <mu, S''[phi]>  for test density phi.
% Use phi = mu itself.
inner_Qmu_mu   = sum(uQ .* (wts_mu.'));   % not really inner product but indicative
inner_Spp_mu_mu = sum((kerspp.eval(chnkrhalf, t)*wts_mu) .* (wts_mu.'));

% boundary evaluation
% <Q[mu], mu>_Gamma  =  <mu, S''[mu]>_Gamma   (self-adjoint in L^2(Gamma))
% Evaluate Q at boundary targets
tb = [];
tb.r = chnkrhalf.r(:,:);
tb.n = chnkrhalf.n(:,:);
tb.d = chnkrhalf.d(:,:);
tb.d2 = chnkrhalf.d2(:,:);

uQ_bdry   = kernq.eval(chnkrhalf, tb) * wts_mu;   % Q[mu] on boundary
uSpp_bdry = kerspp.eval(chnkrhalf, tb) * wts_mu;  % S''[mu] on boundary
bdry_wts  = chnkrhalf.wts(:);

inner_Q   = sum(uQ_bdry   .* (mu .* bdry_wts));
inner_Spp = sum(uSpp_bdry .* (mu .* bdry_wts));

fprintf('[Q vs Spp adjoint]  <Q mu, mu> = %.6e,  <S'''' mu, mu> = %.6e\n', ...
        real(inner_Q), real(inner_Spp));
fprintf('[Q vs Spp adjoint]  rel diff   = %.2e\n', ...
        abs(inner_Q - inner_Spp)/max(abs(inner_Q), abs(inner_Spp)));


%% check Q + Dp
%   (Q + D')[sigma](x)  should equal
%   d_{n_s}^2 u(x) + d_{n_t} d_{n_s} u(x)
%   (Q + D')[sigma]  =  kerndp[sigma]  +  kernq[sigma]
fprintf('\n--- Testing Q + Dp ---\n');

kern_qdp   = kernel('axissymlaplace', 'q_sum_dp',  ndim);
kern_q_sep = kernel('axissymlaplace', 'q',  ndim);

% Evaluate combined kernel
u_qdp = kern_qdp.eval(chnkrhalf, t) * wts;

% Evaluate Q and D' separately and add
u_q_sep  = kern_q_sep.eval(chnkrhalf, t) * wts;
u_dp_sep = kerndp.eval(chnkrhalf, t)     * wts;
u_qdp_sep = u_q_sep + u_dp_sep;

% Compare combined vs sum of parts
err_qdp = log10(abs(u_qdp - u_qdp_sep) ./ max(abs(u_qdp_sep)));
fprintf('[Q+Dp]  max log10 rel-err (combined vs sum-of-parts) = %.2f\n', max(err_qdp));

% Also compare D' part alone against analytic d_{n_t} u
% (FD in target normal direction)
eps_fd = 1e-6;
tp = []; tp.r = t.r + eps_fd*t.n; tp.n = t.n;
tm = []; tm.r = t.r - eps_fd*t.n; tm.n = t.n;

uDp_true_fd = (kernd.eval(chnkrhalf, tp)*wts - kernd.eval(chnkrhalf, tm)*wts)/(2*eps_fd);
err_dp_fd   = log10(abs((u_dp_sep - uDp_true_fd)./max(abs(uDp_true_fd))));
fprintf('[Dp]   max log10 rel-err vs FD = %.2f\n', max(err_dp_fd));

% Cross-check: D'[sigma] vs analytic du/dn_t
dudn_t_true = dudn(t).';   % uses t.n = [0;1]
err_dp_ana  = log10(abs((u_dp_sep - dudn_t_true)./max(abs(dudn_t_true))));
fprintf('[Dp]   max log10 rel-err vs analytic du/dn_t = %.2f\n', max(err_dp_ana));

%% check spp_sum_dp
% S'' + D':  same consistency check
%   (S'' + D')[mu]  =  S''[mu] + D'[mu]
%   d_{n_t}^2 u + d_{n_t} (du/dn_t)  =  d_{n_t}^2 u + d_{n_t} (du/dn_t)
fprintf('\n--- Testing Spp + Dp ---\n');
kern_sppdp    = kernel('axissymlaplace', 'spp_sum_dp', ndim);

% D' with S-density
u_dp_mu    = kerndp.eval(chnkrhalf, t) * wts_mu;
u_spp_mu   = kerspp.eval(chnkrhalf, t) * wts_mu;

% Combined kernel
u_sppdp    = kern_sppdp.eval(chnkrhalf, t) * wts_mu;
u_sppdp_sep = u_spp_mu + u_dp_mu;

err_sppdp   = log10(abs(u_sppdp - u_sppdp_sep) ./ max(abs(u_sppdp_sep)));
fprintf('[Spp+Dp]  max log10 rel-err (combined vs sum-of-parts) = %.2f\n', max(err_sppdp));

% Cross check D' part vs analytic du/dn_t (for S-layer density mu)
err_dp_mu_ana = log10(abs((u_dp_mu - dudn_t_true)./max(abs(dudn_t_true))));
fprintf('[Dp|mu]  max log10 rel-err vs analytic du/dn_t = %.2f\n', max(err_dp_mu_ana));

% Cross check S'' part vs analytic d^2u/dn_t^2
err_spp_ana   = log10(abs((u_spp_mu - d2u_nt2_true)./max(abs(d2u_nt2_true))));
fprintf('[Spp|mu]  max log10 rel-err vs analytic d^2u/dn_t^2 = %.2f\n', max(err_spp_ana));

%% Summary
% plot relative errors for D, S'', D', Q+D', S''+D'
figure(2); clf;
tileplot = tiledlayout(1, 4, 'TileSpacing', 'compact');

titles_list = {'Spp', 'Dp', 'Q+Dp', 'Spp+Dp'};
datas  = {uSpp - d2u_nt2_true, ...
          u_dp_sep - dudn_t_true, ...
          u_qdp - u_qdp_sep, ...
          u_sppdp - u_sppdp_sep};
ref    = {d2u_nt2_true, dudn_t_true, u_qdp_sep, u_sppdp_sep};

for k = 1:4
    ax = nexttile;
    errgrid = nan(size(xx));
    errgrid(in) = log10(abs(datas{k} ./ max(abs(ref{k}))));
    h = pcolor(xx, yy, reshape(errgrid, size(xx)));
    set(h, 'EdgeColor', 'none');
    hold on;
    plot(chnkrhalf); plot(chnkrhalf,'bo'); quiver(chnkrhalf,'r');
    title(titles_list{k});
    clim([-15, -1]);
    colormap(ax, jet(100)); colorbar; axis square;
end

title(tileplot, sprintf('Axisym Laplace kernel tests (n=%d)', ndim));