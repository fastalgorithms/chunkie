function [rhohatinterp,srcinfo,wts] = rhohatInterp(rhohat,rcipsav,ndepth)
%CHNK.RCIP.RHOHATINTERP interpolate the (weight-corrected) coarse level
% density rhohat to the requested depth using the backward recursion
%
% When using an RCIP preconditioner (in the style of eq 34 of the RCIP 
% tutorial), the resulting coarse level density is accurate for
% interactions separated from the vertex *at the coarse scale*. 
%
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

% author: Travis Askham
  
nedge = rcipsav.nedge;
Pbc = rcipsav.Pbc;

starL1 = sort(rcipsav.starL1);
starS = sort(rcipsav.starS);

circL1 = sort(rcipsav.circL1);
circS = sort(rcipsav.circS);

nrho = numel([circS,starS]);
ndens = numel(rhohat)/nrho;
rhohat = reshape(rhohat,[nrho, ndens]);

nsub = rcipsav.nsub;

rhohatinterp = cell(nedge,1);
srcinfo = cell(nedge,1);
wts = cell(nedge,1);

% figure out which edge the indices belong to 
% THIS ASSUMES WE HAVE THE SAME ORDER ON EDGES
circSedge = cell(nedge,1); ncS = numel(circS)/nedge;
circL1edge = cell(nedge,1); ncL1 = numel(circL1)/nedge;
starSedge = cell(nedge,1); nsS = numel(starS)/nedge;
starL1edge = cell(nedge,1); nsL1 = numel(starL1)/nedge;
for j = 1:nedge
    circSedge{j} = circS((j-1)*ncS+1:j*ncS);
    circL1edge{j} = circL1((j-1)*ncL1+1:j*ncL1);
    starSedge{j} = starS((j-1)*nsS+1:j*nsS);
    starL1edge{j} = starL1((j-1)*nsL1+1:j*nsL1);
end

if nargin < 3
    ndepth = nsub;
end

if ndepth > nsub
    msg = "depth requested deeper than RCIP recursion performed\n " + ...
        "going to depth %d instead";
    warning(msg,nsub);
    ndepth = nsub;
end

savedepth = rcipsav.savedepth;

if ndepth <= savedepth
    
    % all relevant quantities are stored, just run backward recursion
    
    rhohat0 = rhohat;
    cl = rcipsav.chnkrlocals{nsub};
    wt = weights(cl);

    for j = 1:nedge
        rhohatinterp{j} = rhohat0(circSedge{j},:);
        srcinfo{j}.r = cl.rstor(:,circL1edge{j});
        srcinfo{j}.d = cl.dstor(:,circL1edge{j});
        srcinfo{j}.d2 = cl.d2stor(:,circL1edge{j});
        srcinfo{j}.n = cl.nstor(:,circL1edge{j});
        wts{j} = wt(circL1edge{j});
    end
    
    R0 = rcipsav.R{nsub+1};
    for i = 1:ndepth
        R1 = rcipsav.R{nsub-i+1};
        MAT = rcipsav.MAT{nsub-i+1};
        rhotemp = R0\rhohat0;
        rhohat0 = R1*(Pbc*rhotemp(starS,:) - MAT*rhohat0(circS,:));
        if i == ndepth
            wt = weights(cl);
            for j = 1:nedge
                rhohatinterp{j} = [rhohatinterp{j}; rhohat0([starSedge{j} circSedge{j}],:)];
                srcinfo{j}.r = [srcinfo{j}.r, cl.rstor(:,starL1edge{j})];
                srcinfo{j}.d = [srcinfo{j}.d, cl.dstor(:,starL1edge{j})];
                srcinfo{j}.d2 = [srcinfo{j}.d2, cl.d2stor(:,starL1edge{j})];
                srcinfo{j}.n = [srcinfo{j}.n, cl.nstor(:,starL1edge{j})];                
                wts{j} = [wts{j}, wt(starL1edge{j})];
            end
        else
            cl = rcipsav.chnkrlocals{nsub-i};
            wt = weights(cl);
            for j = 1:nedge
                rhohatinterp{j} = [rhohatinterp{j}; rhohat0(circSedge{j},:)];
                srcinfo{j}.r = [srcinfo{j}.r, cl.rstor(:,circL1edge{j})];
                srcinfo{j}.d = [srcinfo{j}.d, cl.dstor(:,circL1edge{j})];
                srcinfo{j}.d2 = [srcinfo{j}.d2, cl.d2stor(:,circL1edge{j})];
                srcinfo{j}.n = [srcinfo{j}.n, cl.nstor(:,circL1edge{j})];
                wts{j} = [wts{j}, wt(circL1edge{j})];
            end
        end
        R0 = R1;
    end
    
end
