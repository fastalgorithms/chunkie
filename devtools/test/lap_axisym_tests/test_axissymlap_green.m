% Test values and gradient of Green's third identity for Laplace eq. when 
% the target is on the boundary for the mth mode: 
% 1/2 u = D u - S du/dn
% and
% 1/2 u' = D' u - S' du/dn

% using chunkies hypersingular GGQ because chunkerkerneval only works for
% off/near boundary evaluation

clearvars; 
close all; 
format long e;
%addpaths_loc();

%% geometry
[chnkr] = get_torus_geometry();
charge = [1.0;0.0;-2.0]; % source (axisymm data)

npts = chnkr.npt; % total number of points in discretization
src = chnkr.r(:,:); % generating curve
n_src = chnkr.n(:,:); % normals

% plot geometry
%plot(chnkr);

% convert to cylindrical coordinates
charge_cyl = [sqrt(charge(1)^2 + charge(2)^2);atan2(charge(2),charge(1));charge(3)];

p_modes = 5; % total number of positive fourier modes
n_modes = 2*p_modes + 1; % total number of fourier modes (must be odd for pos/0/neg)
n_angles = n_modes; % number of angles/rotations
m = 3; % mode to test (should be positive and < p_modes, m=1 => 0th mode)

%% compute u and du/dn for densities
u = zeros(n_angles,npts);
dudn = zeros(n_angles,npts);
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

    rvec = pts - charge;
    r = vecnorm(rvec);
    r3 = vecnorm(rvec).^3;
    
    % construct densities
    dudn(i,:) = sum(rvec .* n_pts, 1) ./ (4*pi*r3);
    u(i,:) = 1./(4*pi*r);
end

% Reorder FFT output to match
modes = -p_modes:p_modes;

u_fft = fft(u, n_modes, 1) / n_angles; % FFT (normalized)
dudn_fft = fft(dudn, n_modes, 1) / n_angles; % FFT (normalized)

u_m = fftshift(u_fft, 1);  % puts negative freqs first
dudn_m = fftshift(dudn_fft, 1);  % puts negative freqs first

%% evaluate operators
opts = [];
opts.rcip = false;
opts.forcesmooth = false;
opts.l2scale = false;

Sdudn = 0;
Spdudn = 0;
Du = 0;
Dpu = 0;
origin = [0,0];

% define operators
S = kernel('axissymlaplace','s',m);
Sp = kernel('axissymlaplace','sprime',m);
D = kernel('axissymlaplace','d',m);
Dp = kernel('axissymlaplace','dprime',m);

% evaluate at target
S_m = chunkermat(chnkr, S, opts);
Sp_m = chunkermat(chnkr, Sp, opts);
D_m = chunkermat(chnkr, D, opts);
Dp_m = chunkermat(chnkr, Dp, opts);

% the correct index for the mth mode
i = p_modes + m;

% convolve
Sdudn_m = S_m*dudn_m(i,:)';
Spdudn_m = Sp_m*dudn_m(i,:)';
Du_m = D_m*u_m(i,:)';
Dpu_m = Dp_m*u_m(i,:)';

% check if it satisfies greens identity
val_error = norm(0.5*u_m(i,:)' + Sdudn_m - Du_m);
der_error = norm(0.5*dudn_m(i,:)' + Spdudn_m - Dpu_m);

% display error
disp(['mode:' num2str(modes(i))]);
disp(['green identity error = ' num2str(val_error)]);
disp(['derivative error = ' num2str(der_error)]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [chnkobj] = get_torus_geometry()
    pref = [];
    pref.k = 16; % points per chunk
    %pref.nchmax = 6;

    cparams = [];
    cparams.eps = 1.0e-10;
    %cparams.nover = 1;
    cparams.ifclosed = true;
    cparams.ta = 0;
    cparams.tb = 2*pi;
    cparams.maxchunklen = 2;
    %cparams.nchmin = 4;

    ctr = [3 0];
    narms = 0;
    amp = 0.25;

    chnkobj = chunkerfunc(@(t) starfish(t, narms, amp, ctr), cparams, pref); 
    chnkobj = sort(chnkobj);
end
