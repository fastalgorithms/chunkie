function [rhohatinterp,srcinfo,wts] = rhohatInterp(rhohat,rcipsav,ndepth)
%CHNK.RCIP.RHOHATINTERP interpolate the (weight-corrected) coarse level
% density rhohat to the requested depth using the backward recursion
%
% When using an RCIP preconditioner (in the style of eq 34 of the RCIP 
% tutorial), the resulting coarse level density is 
% accurate for interactions separated from the vertex 
% *at the coarse scale*. 
% By running the recursion for the RCIP preconditioner backwards, an 
% equivalent density can be reconstructed on the fine mesh which is 
% appropriate for closer interactions (at a distance well-separated from
% the vertex *at the finest level in the mesh*).
%
% - Let Gamma_i be the portion of the boundary in the vicinity of
% the vertex at level i (with level nsub the coarse level portion of the
% boundary and level 1 the lowest level). 
% - Let rhohat_i be the weight corrected density on a type c mesh on 
% Gamma_i 
% - Let the matrix I+K_i be a discretization of the integral operator on a
% type b mesh on Gamma_i
% - Let starL denote the "bad indices" of K_i and circL the "good indices";
% the bad correspond to the pairs of small close panels on a type b mesh
% - Let starS denote the "bad indices" of R_i and circS the "good indices";
% the bad correspond to the close panels on a type c mesh
% - Let P be the non-trivial part of an interpolator from type c to type b
% meshes.
% 
% Then, to high precision 
%
% rhohat_{i-1} = R_{i-1}[ [P 0] R_i^(-1) rhohat_i -
%                            (I+K_i)(starL,circL)rhohat_i(circS)]
%
% INPUT
%
% rhohat - entries of density rhohat corresponding to RCIP block in matrix.
%          These come from the two nearest chunks to the corner and 
%          should be ordered so that rhohat(starS) and rhohat(circS)
%          correspond to the density at the bad and good points (near and
%          far chunk) in the correct order (accounting for chunk
%          orientation)
% rcipsav - the quantities obtained from a previous call to Rcompchunk 
% ndepth - (default: rcip.nsub) compute up to rhohat_{nsub-ndepth+1} 
%
% OUTPUT
%
% rhohatinterp - array of interpolated values, 
% [rhohat_nsub(circS); rhohat_{nsub-1}(circS); rhohat_{nsub-2}(circS) ... 
%                      rhohat_1(circS); rhohat_1(starS)]
%                note that these are mostly the "good" indices except at
%                the lowest level (naturally)
% srcinfo - struct, containing points (srcinfo.r), normals (srcinfo.n), etc
%                 in the same order as rhohatinterp.
% wts - a set of weights for integrating functions sampled at these
%           points


rhohatinterp = [];
srcinfo = [];
wts = [];

k = rcipsav.k;
ndim = rcipsav.ndim;
nedge = rcipsav.nedge;
Pbc = rcipsav.Pbc;
PWbc = rcipsav.PWbc;
starL = rcipsav.starL;
starL1 = rcipsav.starL1;
starS = rcipsav.starS;
circL = rcipsav.circL;
circL1 = rcipsav.circL1;
circS = rcipsav.circS;
ilist = rcipsav.ilist;
nsub = rcipsav.nsub;

if nargin < 3
    ndepth = nsub;
end

if ndepth > nsub
    msg = "depth requested deeper than RCIP recursion performed\n " + ...
        "going to depth %d instead";
    warning(msg,nsub);
    ndepth = nsub;
end

savedeep = rcipsav.savedeep;

if savedeep
    
    % all relevant quantities are stored, just run backward recursion
    
    rhohat0 = rhohat;
    rhohatinterp = [rhohatinterp; rhohat0(circS)];
    r = [];
    d = [];
    d2 = [];
    n = [];
    h = [];
    cl = rcipsav.chnkrlocals{nsub};
    wt = weights(cl);
    r = [r, cl.rstor(:,circL1)];
    d = [d, cl.dstor(:,circL1)];
    d2 = [d2, cl.d2stor(:,circL1)];
    n = [n, cl.nstor(:,circL1)];
    wts = [wts; wt(circL1(:))];
    
    R0 = rcipsav.R{nsub+1};
    for i = 1:ndepth
        R1 = rcipsav.R{nsub-i+1};
        MAT = rcipsav.MAT{nsub-i+1};
        rhotemp = R0\rhohat0;
        rhohat0 = R1*(Pbc*rhotemp(starS) - MAT*rhohat0(circS));
        if i == ndepth
            rhohatinterp = [rhohatinterp; rhohat0];
            r = [r, cl.rstor(:,starL1)];
            d = [d, cl.dstor(:,starL1)];
            d2 = [d2, cl.d2stor(:,starL1)];
            n = [n, cl.nstor(:,starL1)];
            wt = weights(cl);

            wts = [wts; wt(starL1(:))];
        else
            cl = rcipsav.chnkrlocals{nsub-i};
            rhohatinterp = [rhohatinterp; rhohat0(circS)];
            r = [r, cl.rstor(:,circL1)];
            d = [d, cl.dstor(:,circL1)];
            d2 = [d2, cl.d2stor(:,circL1)];
            n = [n, cl.nstor(:,circL1)];
            wt = weights(cl);
            wts = [wts; wt(circL1(:))];
        end
        R0 = R1;
    end
    
    srcinfo.r = r;
    srcinfo.d = d;
    srcinfo.d2 = d2;
    srcinfo.n = n;

end
