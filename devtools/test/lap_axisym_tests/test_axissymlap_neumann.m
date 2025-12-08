%{
- 3D Laplace's equation
- Neumann boundary condition
- Torus boundary
- Single and Double layer potential representation
- mth mode (non-axisymmetric B.C.)
%}

clearvars; 
close all;
format long e;

%% geometry
% target is where we evaluate the solution
[chnkr,target,charge1,charge2] = get_torus_geometry();

npts = chnkr.npt; % total number of points in discretization
src = chnkr.r(:,:); % generating curve
n_src = chnkr.n(:,:); % normals

% plot geometry
%plot(chnkr);

p_modes = 3; % number of positive fourier modes
n_modes = 2*p_modes + 1; % number of fourier modes (must be odd for pos/0/neg)
n_angles = n_modes; % number of angles/rotations
strength = 1.0; % strength of charges

%% discretization
% compute f (boundary condition)
f = zeros(n_angles,npts);
for i=1:n_angles
    % get polar/cartesian coordinates
    theta = (i-1)*2*pi/n_angles;
    
    x=src(1,:).*cos(theta);
    y=src(1,:).*sin(theta);
    z=src(2,:);
    pts=[x;y;z];

    n_x=n_src(1,:).*cos(theta);
    n_y=n_src(1,:).*sin(theta);
    n_z=n_src(2,:);
    n_pts = [n_x;n_y;n_z];

    % set up neumann boundary data (with both charges)
    rvec1 = pts - charge1;
    r3_1 = vecnorm(rvec1).^3;
    fn1 = sum(rvec1 .* n_pts, 1) ./ (4*pi*r3_1);

    rvec2 = pts - charge2;
    r3_2 = vecnorm(rvec2).^3;
    fn2 = -sum(rvec2 .* n_pts, 1) ./ (4*pi*r3_2);

    f(i,:) = strength*(fn1 + fn2);
end

% Reorder FFT output to match
modes = -p_modes:p_modes;
f_fft = fft(f, n_modes, 1) / n_angles; % FFT (normalized)
f_m = fftshift(f_fft, 1);  % puts negative freqs first

%% solve
% solve the integral equation for each fourier mode
opts = [];
opts.rcip = false;
opts.forcesmooth = false;
opts.l2scale = false;

sigma1_m = zeros(n_modes,npts); % single layer density
sigma2_m = zeros(n_modes,npts); % double layer density
origin = [0,0];
for i=1:n_modes
    m = abs(modes(i))+1;
    Sp = kernel('axissymlap','sprime',m);
    Dp = kernel('axissymlap','dprime',m);

    % Build the system matrix
    Sp_m = chunkermat(chnkr, Sp, opts) - 0.5*eye(npts);
    Dp_m = chunkermat(chnkr, Dp, opts);

    % Enforce zero-mean constraint for compatability condition
    Sp_m = Sp_m + onesmat(chnkr);
    Dp_m = Dp_m + onesmat(chnkr);

    % Solve the linear system
    rhs = f_m(i,:)';
    sigma1_m(i,:) = gmres(Sp_m, rhs, [], 1e-12, npts);
    sigma2_m(i,:) = gmres(Dp_m, rhs, [], 1e-12, npts);
end

%% solution building
opts.verb = false;
opts.quadkgparams = {'RelTol', 1e-10, 'AbsTol', 1.0e-10};

% target in cylindrical coordinates (r,theta,z)
target_cyl = [sqrt(target(1)^2 + target(2)^2);atan2(target(2),target(1));target(3)];
target_new = [target_cyl(1);target_cyl(3)];

u1_sol = 0; % single layer solution
u2_sol = 0; % double layer solution
for i=1:n_modes
    m = abs(modes(i))+1;
    S = kernel('axissymlap','s',m);
    D = kernel('axissymlap','d',m);
    
    % evaluation of operators
    u1_m = chunkerkerneval(chnkr, S, sigma1_m(i,:), target_new, opts);
    u2_m = chunkerkerneval(chnkr, D, sigma2_m(i,:), target_new, opts);

    % fourier composition
    u1_sol = u1_sol + real(u1_m * exp(1i * modes(i) * target_cyl(2)));
    u2_sol = u2_sol + real(u2_m * exp(1i * modes(i) * target_cyl(2)));
end

% compute the exact solution explicitly
% only unique up to a constant since neumann BC
r1 = norm(target - charge1);
r2 = norm(target - charge2);
u_true = strength*1.0/(4*pi)*(1/r1 - 1/r2);

% compute the error of 1st/2nd kind integral equation
err1 = norm(u1_sol-u_true)
err2 = norm(u2_sol-u_true)




%% gemotry functions
function [chnkobj,target,charge1,charge2] = get_torus_geometry()
    pref = [];
    pref.k = 16; % points per chunk
    %pref.nchmax = 2;

    cparams = [];
    %cparams.eps = 1.0e-10;
    %cparams.nover = 1;
    cparams.ifclosed = true;
    cparams.ta = 0;
    cparams.tb = 2*pi;
    cparams.maxchunklen = 2;
    %cparams.nchmin = 8;

    ctr = [3 0];
    narms = 0;
    amp = 0.25;

    chnkobj = chunkerfunc(@(t) starfish(t, narms, amp, ctr), cparams, pref); 
    chnkobj = sort(chnkobj);

    target = [3;0.0;-0.7];
    charge1 = [1.0;0.0;3.0];
    charge2 = [1.0;0.0;-3.0];
end
