clear 
clc

kvec = [1;-1.5];

%

zk = sqrt(norm(kvec));               % our k (wave number)
nu = 1/3;
cparams = [];

cparams.eps = 1e-6;
%cparams.nover = 0;
cparams.maxchunklen = 4./zk;       % setting a chunk length helps when the
                                    % frequency is known'



chnkr = chunkerfunc(@(t) starfish(t), cparams);

figure(1)                                                   % plot the chunker-object (supposed to be a circle centered at 1 with radius 1)
clf
plot(chnkr, '-x')
title('Chunkr object')
hold on
quiver(chnkr)
axis equal
drawnow()

coefs = [nu; 0];
opts = [];
opts.sing = 'log';
fkern =  @(s,t) chnk.helm2d.kern(zk, s, t, 'free plate first part', coefs);        % build the desired kernel


fkern1 = @(s,t) chnk.helm2d.kern(zk, s, t, 'free plate hilbert subtract', coefs);                   % hilbert subtraction kernels in K11
fkern2 = @(s,t) chnk.helm2d.kern(zk, s, t, 'free plate coupled hilbert', coefs);   
hilbert = @(s,t) chnk.lap2d.kern(s, t, 'hilb');
double = @(s,t) chnk.lap2d.kern(s,t, 'd');

fkern3 = @(s,t) chnk.helm2d.kern(zk, s, t, 'free plate K21 first part', coefs);                     % singularity subtration kernel in K21 (including swapping its Asmyptotics expansions)

fkern4 = @(s,t) chnk.helm2d.kern(zk, s, t, 'free plate K21 second part', coefs);                    % kernels in K21 needs to multiply by curvature

fkern5 = @(s,t) chnk.helm2d.kern(zk, s, t, 'free plate K21 hilbert part', coefs);                   % kernels in K21 coupled with hilbert transforms and needs to multiply by curvature

fkern6 = @(s,t) chnk.helm2d.kern(zk, s, t, 'free plate K22 second part', coefs);                    % kernels in K22 needs to multiply by curvature


start = tic;
sysmat = chunkermat(chnkr,fkern, opts);
sysmat1 = chunkermat(chnkr, fkern1, opts);
sysmat2 = chunkermat(chnkr, fkern2, opts);
K21 = chunkermat(chnkr, fkern3, opts);
K21second = chunkermat(chnkr, fkern4, opts);
K21hilbert = chunkermat(chnkr, fkern5, opts);

K22second = chunkermat(chnkr, fkern6, opts);

D = chunkermat(chnkr, double, opts);
t1 = toc(start);
fprintf('%5.2e s : time to assemble matrix\n',t1)

opts2 = [];
opts2.sing = 'pv';

start = tic;
H = chunkermat(chnkr, hilbert, opts2);                                              % Assemble hilbert transforms
t1 = toc(start);
fprintf('%5.2e s : time to assemble matrix\n',t1)


hilb = sysmat1*H - ((1+nu)/2).*(D*D)- ((1+nu)*nu/2).*(D*D);
hilb2 = sysmat2*H + K21hilbert*H;

mat1 =  sysmat(1:2:end, 1:2:end);
mat4 =  sysmat(2:2:end, 2:2:end);



sysmat(1:2:end, 1:2:end) = mat1 + hilb;
sysmat(2:2:end, 1:2:end) = K21 +  hilb2 + K21second;
sysmat(2:2:end, 2:2:end) = mat4 + K22second;




A = [-1/2 + (1/8)*(1+nu).^2, 0; 0, 1/2];                                     % jump matrix (for interior problem)

M = kron(eye(chnkr.npt), A);

lhs =  M + sysmat ;

[hess, third, ~] = chnk.helm2d.thirdforth_derivatives(zk, [0;0], chnkr.r);
[hessK, thirdK, ~] = chnk.helm2d.thirdforth_derivatives(zk*(1i), [0;0], chnkr.r);

nx = chnkr.n(1,:).'; 
ny = chnkr.n(2,:).';

dx = chnkr.d(1,:).';
dy = chnkr.d(2,:).';

ds = sqrt(dx.*dx+dy.*dy);
taux = (dx./ds);                                                                       % normalization
tauy = (dy./ds);

firstbc = 1/(2*zk^2).*(hess(:, :, 1).*(nx.*nx) + hess(:, :, 2).*(2*nx.*ny) + hess(:, :, 3).*(ny.*ny))-...
           1/(2*zk^2).*(hessK(:, :, 1).*(nx.*nx) + hessK(:, :, 2).*(2*nx.*ny) + hessK(:, :, 3).*(ny.*ny))+...
           coefs(1)/(2*zk^2).*(hess(:, :, 1).*(taux.*taux) + hess(:, :, 2).*(2*taux.*tauy) + hess(:, :, 3).*(tauy.*tauy))-...
           coefs(1)/(2*zk^2).*(hessK(:, :, 1).*(taux.*taux) + hessK(:, :, 2).*(2*taux.*tauy) + ...
           hessK(:, :, 3).*(tauy.*tauy));

secondbc = -1./(2*zk^2).*(third(:, :, 1).*(nx.*nx.*nx) + third(:, :, 2).*(3*nx.*nx.*ny) +...
       third(:, :, 4).*(3*nx.*ny.*ny) + third(:, :, 6).*(ny.*ny.*ny)) + ...
        1./(2*zk^2).*(thirdK(:, :, 1).*(nx.*nx.*nx) + thirdK(:, :, 2).*(3*nx.*nx.*ny)+...
        thirdK(:, :, 4).*(3*nx.*ny.*ny) + thirdK(:, :, 6).*(ny.*ny.*ny)) -...
        (2-coefs(1))/(2*zk^2).*(third(:, :, 1).*(taux.*taux.*nx) + third(:, :, 2).*(taux.*taux.*ny + 2*taux.*tauy.*nx) +...
        third(:, :, 4).*(2*taux.*tauy.*ny+ tauy.*tauy.*nx) +...
        + third(:, :, 6).*(tauy.*tauy.*ny)) + ...
        (2-coefs(1))/(2*zk^2).*(thirdK(:, :, 1).*(taux.*taux.*nx) + thirdK(:, :, 2).*(taux.*taux.*ny + 2*taux.*tauy.*nx) +...
        thirdK(:, :, 4).*(2*taux.*tauy.*ny+ tauy.*tauy.*nx) +...
        + thirdK(:, :, 6).*(tauy.*tauy.*ny)) +...
        (1-coefs(1)).*(1/(2*zk^2).*(hess(:, :, 1).*taux.*taux + hess(:, :, 2).*(2*taux.*tauy) + hess(:, :, 3).*tauy.*tauy)-...
         1/(2*zk^2).*(hessK(:, :, 1).*taux.*taux + hessK(:, :, 2).*(2*taux.*tauy) + ...
        hessK(:, :, 3).*tauy.*tauy)-...
        (1/(2*zk^2).*(hess(:, :, 1).*nx.*nx + hess(:, :, 2).*(2*nx.*ny) + hess(:, :, 3).*ny.*ny)-...
       1/(2*zk^2).*(hessK(:, :, 1).*nx.*nx + hessK(:, :, 2).*(2*nx.*ny) + hessK(:, :, 3).*ny.*ny)));

[nt, ~] = size(sysmat);

rhs = zeros(nt, 1); 
rhs(1:2:end) = firstbc ; 
rhs(2:2:end) = secondbc;

tic
%sol = gmres(lhs, rhs, [], 1e-13, 200);
sol = lhs\rhs;
toc;

rho1 = sol(1:2:end);                                    % first density
rho2 = sol(2:2:end);        


xs = -3:0.01:3;                                    % generate some targets
ys = -3:0.01:3;
[X,Y] = meshgrid(xs, ys);
targets = [X(:).'; Y(:).'];
[~,na] = size(targets);

tic
in = chunkerinterior(chnkr, targets); 
out = ~in;
trad = targets(1,:).^2+targets(2,:).^2;
%out = out.*(sqrt(trad.')>2);
out = find(out);
toc


ikern1 = @(s,t) chnk.helm2d.kern(zk, s, t, 'free plate eval first', coefs);                              % build the kernel of evaluation          
ikern2 = @(s,t) chnk.helm2d.kern(zk, s, t, 'free plate eval second');
ikern3 = @(s,t) chnk.helm2d.kern(zk, s, t, 'free plate eval first hilbert',coefs);

opts_off = [];
opts_off.forcesmooth = true;
opts_off.accel = false;
%coupled = chunkerkernevalmat(chnkr, ikern3, targets(:, out), opts_off);

start1 = tic;
Hrho = H*rho1;
Dsol1 = chunkerkerneval(chnkr, ikern1,rho1, targets(:, out),opts_off);
Dsol2 = chunkerkerneval(chnkr, ikern3, Hrho, targets(:, out),opts_off);
Dsol3 = chunkerkerneval(chnkr, ikern2, rho2, targets(:,out),opts_off);
Dsol = Dsol1 + Dsol2 + Dsol3;
t2 = toc(start1);
fprintf('%5.2e s : time for kernel eval (for plotting)\n',t2)

true_sol = zeros(na, 1);
utarg = zeros(na, 1);

[val,~] = chnk.helm2d.green(zk,[0;0],targets(:,out));        % Hankel part

zkimag = (1i)*zk;
[valK,~] = chnk.helm2d.green(zkimag,[0;0], targets(:,out));    % modified bessel K part

trueval = 1/(2*zk^2).*val - 1/(2*zk^2).*valK;


utarg(out) = Dsol;
true_sol(out) = trueval;


uerr = utarg - true_sol;
uerr = reshape(uerr,size(X));
figure(2)
h = pcolor(X,Y,log10(abs(uerr)));
set(h,'EdgeColor','None'); hold on;
title("Absolute error (free plate kernel on a circle)", 'FontSize',16)
plot(chnkr,'w-','LineWidth',2);

colorbar

