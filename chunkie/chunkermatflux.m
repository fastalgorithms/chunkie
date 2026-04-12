function fluxes = chunkermatflux(chnkobj,kern,dens,fun,rcipsav,opts)
%%%CHUNKERMATFLUX evaluate integral of layer potential on the boundary.
% This is particularly useful when RCIP is used, which prevents the naive
% use of chunkermat
% 
% Syntax: 
%   fluxes = chunkermat(chnkr,kerns,sol)
%   fluxes = chunkermat(chnkr,kerns,sol,fun)
%   fluxes = chunkermat(chnkr,kerns,sol,fun,rcipsav)
% Input:
%   chnkobj - chunker object describing boundary
%   kern - kernel function or matrix of kernels. 
%   dens - layer potential density 
%   fun - if provided, this will routine will compute int_gamma
%       fun*K[dens]. fun will be interpolated down to the fine grid. 
%       Default is fun = 1
%   rcipsav - precomputed structure of rcip data at corners (see
%       chunkermat)
%   opts  - options structure. available options (default settings)
%       opts.tot - (false) if true, return the sum over all edges as a
%                  number, rather than individual edges
%       opts.matopts - opts struct to pass to global chunkermat
% Output:
%   fluxes - integral of fun(x,y)*K[dens] over each edge of the chunker
%
% TODO: what if the boundary is just a chunker? What about if it is an
% array of chunkers?
% Vector valued kernels/densities
% What if RCIP was off?
% Test computing int_omega u.
% Test with a chunkgraph with 3 incident edges.
% Compute RCIP when we have to apply a singular op to a smooth rhs, e.g.
%   Sprime on a polygon applied to piecewise smooth density
% pass opts to chunkermat, pass ilist to chunkermat?
% change chunkermat to chunkermatapply?
% return local chunkgraphs, local densities, and local layer potentials

if (class(chnkobj) == "chunker")
    chnkobj = tochunkgraph(chnkobj);
elseif~(class(chnkobj) == "chunkgraph")
    msg = "CHUNKERMATFLUX: first input is not a chunker or chunkgraph object";
    error(msg)
end

npts = [chnkobj.echnks.npt];
nedge = length(npts);
nvert = size(chnkobj.verts,2);

if ~isa(kern,'kernel')
    try 
        kern = kernel(kern);
    catch
        error('CHUNKERMAT: second input kern not of supported type');
    end
end

if nargin < 4
    fun = ones(1,chnkobj.npt);
end
if isempty(fun)
    fun = ones(1,chnkobj.npt);
end

if nargin < 5
    rcipsav = [];
end

matopts = [];
if isfield(opts,'matopts')
    matopts = opts.matopts;
end

matopts.rcip = false;
sysmat0 = chunkermat(chnkobj,kern,matopts);


% get starting index for each edge
idstart = [1,cumsum(npts)+1];

ids_rcip = cell(2,nedge);
ids_ignore = cell(2,nedge);

% create list of interactions done incorrectly by chunkermat
for i = 1:nedge
    ids_rcip{1,i} = idstart(i)+(1:2*16)-1;
    ids_rcip{2,i} = flip(idstart(i+1)-(1:2*16));

    ids_ignore{1,i} = idstart(i)+(1:3*16)-1;
    ids_ignore{2,i} = flip(idstart(i+1)-(1:3*16));
end

% zero out interactions done incorrectly by chunkermat
for i = 1:nvert
    vstruc = chnkobj.vstruc{i};
    for j = 1:length(vstruc{1})
        jedge = vstruc{1}(j);
        for k = 1:length(vstruc{1})
        kedge = vstruc{1}(k);
        sysmat0(ids_ignore{1.5+vstruc{2}(j)/2,jedge}, ids_rcip{1.5+vstruc{2}(k)/2,kedge}) = 0;
        end
    end
end

% compute the contribution to the fluxes done correctly by chunkermat
val = sysmat0*dens;

ids = [];
fluxes = zeros(1,nedge);
for i = 1:nedge
    ids = edgeids(chnkobj,i);
    fluxes(i) = sum(val(ids).*chnkobj.wts(ids).');
end

% now go back and add the vertex contributions
nverts = size(chnkobj.verts,2);

% TODO: grab this from rcipsav
ndepth = 40;
for ivert = 1:nverts
    starindtmp = rcipsav{ivert}.starind;
    
    solhat = dens(starindtmp);
    [solhatinterpcell,srcinfocell,wtscell] = chnk.rcip.rhohatInterp(solhat,rcipsav{ivert},ndepth);

    % construct panels and density for this vertex
    solinterp = [];
    % solhatinterp2 = [];
    chnkrs = cell(1,length(srcinfocell));
    for j = 1:length(srcinfocell)
        solhat = []; solhat.val = solhatinterpcell{j}(:).';

        chnkloc = chnkobj.echnks(chnkobj.vstruc{ivert}{1}(j));
        if chnkobj.vstruc{ivert}{2}(j) == -1
            srcinfocell{j} = flip_struct(srcinfocell{j});
            solhat = flip_struct(solhat);

            rloc = flip(chnkloc.r(:,:,3),2) - rcipsav{ivert}.ctr(:,j);
            dloc = flip(chnkloc.d(:,:,3),2);
            d2loc = flip(chnkloc.d2(:,:,3),2);
        else
            rloc = chnkloc.r(:,:,end-2) - rcipsav{ivert}.ctr(:,j);
            dloc = chnkloc.d(:,:,end-2);
            d2loc = chnkloc.d2(:,:,end-2);
        end
        
        src = [];
        src.r = reshape([rloc,srcinfocell{j}.r],2,16,[]);
        src.d = reshape([dloc,srcinfocell{j}.d],2,16,[]);
        src.d2 = reshape([d2loc,srcinfocell{j}.d2],2,16,[]);
        chnkrj = chunkerpoints(src);
        chnkrs{j} = chnkrj;

        % solinterp = [solinterp; solhat.val(:)];
        solinterp = [solinterp; zeros(16,1);solhat.val(:)];
    end

    % construct local chunkgraph for this edge (centered at [0;0])
    verts = [[0;0]];
    edge2verts = [[NaN;NaN],[NaN;NaN]];
    cgrphloc = chunkgraph(verts,edge2verts,chnkrs);

    % get local kernels
    clear fkernloc
    if (size(kern) == 1)
        fkernloc = kern;
    else
        fkernloc(length(srcinfocell),length(srcinfocell)) = kernel();
        fkernloc(:,:) = kern(chnkobj.vstruc{ivert}{1},chnkobj.vstruc{ivert}{1});
    end 

    % get local chunkermat
    opts = []; opts.rcip = false; %opts.adaptive_correction = true;
    sysmatloc = chunkermat(cgrphloc,fkernloc,opts);

    val_loc = sysmatloc * solinterp;

    % val_loc = val_loc + solinterp;

    % compute local contribution to the fluxes
    for j = 1:length(srcinfocell)
        ids = edgeids(cgrphloc,j);
        flux_loc = sum(val_loc(ids).*cgrphloc.wts(ids).');

        fluxes(chnkobj.vstruc{ivert}{1}(j)) = ... 
            fluxes(chnkobj.vstruc{ivert}{1}(j)) + flux_loc;
    end
end

% add up all fluxes, if requested.
if isfield(opts,'tot')
    if opts.tot
        fluxes = sum(fluxes);
    end
end
end