function [proxy_dens,bragg_coef,interface_dens] = solve_sys(full_sys,rhs,KK)

nwall=120; % number of points in each wall segment 
nproxy = 160;

% Blocksolve the system.
invA = inv(full_sys.Amat);
% 
AinvB = invA*full_sys.Bmat;

Schur_top=full_sys.Dmat-full_sys.Cmat*AinvB;
Schur_bot=full_sys.Vmat-full_sys.Zmat*AinvB;

num_layer = 2;
Schur3=[Schur_top, zeros(num_layer*nwall*2,2*KK);...
    Schur_bot, full_sys.Wmat];

[SU,SS,SV] = svd(Schur3);

ind = find(diag(SS)>1e-13);
dia = diag(SS);
dia = dia(ind);

SU = SU(:,ind);
SS = diag(1./dia);
SV = SV(:,ind);

% A^-1*f
Ainvf = invA*rhs;

rhs_schur3=-[full_sys.Cmat*Ainvf;...
    full_sys.Zmat*Ainvf];

% proxy circle charges
coeff=SV*(SS*(SU'*(rhs_schur3)));

% Rayleigh Bragg coefficients (top and bottom)
bragg_coef=coeff(num_layer*nproxy+1:num_layer*nproxy+2*KK);
% proxy circle charges
proxy_dens=coeff(1:num_layer*nproxy);
% Interface densities
interface_dens=Ainvf-AinvB*proxy_dens;




return