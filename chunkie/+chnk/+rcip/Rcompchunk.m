function [R,rcipsav]=Rcompchunk(chnkr,iedgechunks,fkern,ndim, ...
    Pbc,PWbc,nsub,starL,circL,starS,circS,ilist,starL1,circL1,...
    sbclmat,sbcrmat,lvmat,rvmat,u,opts)
%CHNK.RCIP.Rcompchunk carry out the forward recursion for computing
% the preconditioner R where geometry is described as a chunker
%
% This routine is not intended to be user-callable 
%
% Function is passed as a handle, number of equations is given by
% ndim
%
% Kernel on input takes in arguments (chnkrlocal,ilistl);
% 
% Note that matrix must be scaled to have identity on the diagonal,
% will not work with scaled version of identity
%

% author: Shidong Jiang
% modified: Jeremy Hoskins, Manas Rachh


k = chnkr.k;  
dim = chnkr.dim;
nedge = size(iedgechunks,2);

glxs = chnkr.tstor;
glws = chnkr.wstor;

% return what's needed to interpolate from coarse

rcipsav = [];
rcipsav.k = k;
rcipsav.ndim = ndim;
rcipsav.nedge = nedge;
rcipsav.Pbc = Pbc;
rcipsav.PWbc = PWbc;
rcipsav.starL = starL;
rcipsav.starL1 = starL1;
rcipsav.starS = starS;
rcipsav.circL = circL;
rcipsav.circL1 = circL1;
rcipsav.circS = circS;
rcipsav.ilist = ilist;
rcipsav.nsub = nsub;

if (nargin < 15 || isempty(sbclmat) || isempty(sbcrmat) || ...
        isempty(lvmat) || isempty(rvmat) || isempty(u))
    [sbclmat,sbcrmat,lvmat,rvmat,u] = chnk.rcip.shiftedlegbasismats(k); 
end

if nargin < 20
    opts = [];
end

savedeep = false;
if isfield(opts,'savedeep')
    savedeep = opts.savedeep;
end

rcipsav.savedeep = savedeep;

if savedeep
    rcipsav.R = cell(nsub+1,1);
    rcipsav.MAT = cell(nsub,1);
    rcipsav.chnkrlocals = cell(nsub,1);
end


% grab only those kernels relevant to this vertex

if(size(fkern)==1)
    fkernlocal = fkern;
else
  fkernlocal(nedge,nedge) = kernel();
    for i=1:nedge
        ici = iedgechunks(1,i);
        for j=1:nedge
            icj = iedgechunks(1,j);
            fkernlocal(i,j) = fkern(ici,icj);
        end
    end

end

rcipsav.fkernlocal = fkernlocal;

% get coefficients of recentered edge chunks and figure out orientation

km1 = k-1;
rcs = zeros(km1,dim,nedge);
dcs = zeros(k,dim,nedge);
d2cs = zeros(k,dim,nedge);
dscal = zeros(nedge,1);
d2scal = zeros(nedge,1);
ctr = zeros(dim,nedge);

ileftright = zeros(nedge,1);
nextchunk = zeros(nedge,1);

for i = 1:nedge
    ic = iedgechunks(1,i);
    ie = iedgechunks(2,i);
    chnkri = chnkr(ic);
    r = chnkri.r(:,:,ie);
    d = chnkri.d(:,:,ie);    
    d2 = chnkri.d2(:,:,ie);    
    il = chnkri.adj(1,ie);
    ir = chnkri.adj(2,ie);
    h = chnkri.h(ie);
    if (il > 0 && ir < 0)
        nextchunk(i) = il;
        ileftright(i) = 1;

        rr = rvmat*(r.'); r = r - rr(:);
        ctr(:,i) = rr;
        rcs(:,:,i) = sbcrmat*(r.');
        dcs(:,:,i) = u*(d.');
        d2cs(:,:,i) = u*(d2.');
        dscal(i) = h*2;
        d2scal(i) = h^2*4; 
    elseif (il < 0 && ir > 0)
        nextchunk(i) = ir;
        ileftright(i) = -1;
        rl = lvmat*(r.'); r = r - rl(:);
        ctr(:,i) = rl;
        rcs(:,:,i) = sbclmat*(r.');
        dcs(:,:,i) = u*(d.');
        d2cs(:,:,i) = u*(d2.');
        dscal(i) = h*2;
        d2scal(i) = h^2*4; 
    else
        error('RCIP: edge chunk not adjacent to one vertex and one neighbor')
    end
end

rcipsav.ctr = ctr;
rcipsav.rcs = rcs;
rcipsav.dcs = dcs;
rcipsav.d2cs = d2cs;
rcipsav.dscal = dscal;
rcipsav.d2scal = d2scal;
rcipsav.ileftright = ileftright;
rcipsav.glxs = glxs;
rcipsav.glws = glws;

pref = []; 
pref.k = k;
pref.nchstor = 5;

R = [];

% size of the system matrix
nsys = 3*k*nedge*ndim;
  
% size of the preconditioner R
nR = 2*k*nedge*ndim;

ts = cell(nedge,1);
chnkrlocal(1,nedge) = chunker();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% begin recursion proper

h0=ones(nedge,1);
for level=1:nsub
    h = h0/2^(nsub-level);

    for i=1:nedge
        if ileftright(i) == -1
            if level == nsub 
                ts{i} =  [0, 0.5, 1]*h(i);
            else
                ts{i} =  [0, 0.5, 1, 2]*h(i);
            end
        else
            if level == nsub
                ts{i} = -[1, 0.5, 0]*h(i);
            else
                ts{i} = -[2, 1, 0.5, 0]*h(i);
            end
            
        end
    end
    % construct local chunks around the corner
    for i=1:nedge
        chnkrlocal(i) = chnk.rcip.chunkerfunclocal(@(t) shiftedcurve(t,rcs(:,:,i),dcs(:,:,i), ...
            dscal(i),d2cs(:,:,i),d2scal(i),ileftright(i)),ts{i},pref,glxs,glws);
    end
    
    % at the top level, append/prepend the next chunk
    if level == nsub
        for i = 1:nedge
            ic = iedgechunks(1,i);
            nc = nextchunk(i);
            nchi = chnkrlocal(i).nch;
            chnkrlocal(i) = chnkrlocal(i).addchunk(1);
            chnkrlocal(i).r(:,:,nchi+1) = chnkr(ic).r(:,:,nc)-ctr(:,i);
            chnkrlocal(i).d(:,:,nchi+1) = chnkr(ic).d(:,:,nc);                
            chnkrlocal(i).d2(:,:,nchi+1) = chnkr(ic).d2(:,:,nc); 
            chnkrlocal(i).h(nchi+1) = chnkr(ic).h(nc);
            if ileftright(i) == -1
                chnkrlocal(i).adj(1,nchi+1) = nchi;
                chnkrlocal(i).adj(2,nchi+1) = -1;
                chnkrlocal(i).adj(2,nchi) = nchi+1;
            else
                chnkrlocal(i).adj(1,nchi+1) = -1;
                chnkrlocal(i).adj(2,nchi+1) = 1;
                chnkrlocal(i).adj(1,1) = nchi+1;
                chnkrlocal(i) = chnkrlocal(i).sort();
            end
            chnkrlocal(i).n = normals(chnkrlocal(i));
            chnkrlocal(i).wts = weights(chnkrlocal(i));
        end
    end
    
% construct the system matrix for local chunks
    if level == 1
        ilistl = [];
    else
        ilistl = ilist;
    end

    % test for opdims ~= [1,1]
    [MAT,opts] = chunkermat(chnkrlocal,fkernlocal,opts,ilistl);
    

%
    MAT = eye(nsys) + MAT;
    if level==1    %  Dumb and lazy initializer for R, for now
  %R=eye(nR); 
        R = inv(MAT(starL,starL));
        if savedeep
            rcipsav.R{1} = R;
        end
    end
    R=chnk.rcip.SchurBana(Pbc,PWbc,MAT,R,starL,circL,starS,circS);   
    if savedeep
        rcipsav.R{level+1} = R;
        rcipsav.MAT{level} = MAT(starL,circL);
        rcipsav.chnkrlocals{level} = merge(chnkrlocal);
    end
end

end

function [r,d,d2] = shiftedcurve(t,rc,dc,scald,d2c,scald2,ilr)
n = length(dc);
nm1 = n-1;
if (ilr == -1)
    tt = 2*t-1;
else
    tt = 2*t+1;
end
pols = lege.pols(tt,nm1); pols = pols.';
r = (t.*(pols(:,1:nm1)*rc)).';
d = (scald*(pols*dc)).';
d2 = (scald2*(pols*d2c)).';

end

