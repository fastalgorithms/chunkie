function fluxes = chunkermatflux(chunkobj,kerns,dens,fun,rcipsav,opts)
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
%   kerns - kernel function or matrix of kernels. 
%   dens - layer potential density 
%   fun - if provided, this will routine will compute int_gamma
%       fun(x,y)*K[dens]. Default is fun(x,y) = 1
%   rcipsav - precomputed structure of rcip data at corners (see
%       chunkermat)
%   opts  - options structure. available options (default settings)
%       opts.tot - (false) if true, return the sum over all edges as a
%                  number, rather than individual edges
% Output:
%   fluxes - integral of fun(x,y)*K[dens] over each edge of the chunker
%
% TODO: what if the boundary is just a chunker? Vector valued
% kernels/densities
% test with fun to give int_omega u.

opts.rcip = false;
sysmat0 = chunkermat(chunkobj,kerns,opts);

npts = [chunkobj.echnks.npt];
nedge = length(npts);
nvert = size(chunkobj.verts,2);

% get starting index for each edge
idstart = [1,cumsum(npts)+1];

ids_rcip = cell(2,nedge);
ids_ignore = cell(2,nedge);

for i = 1:nedge
    ids_rcip{1,i} = idstart(i)+(1:2*16)-1;
    ids_rcip{2,i} = flip(idstart(i+1)-(1:2*16));

    ids_ignore{1,i} = idstart(i)+(1:3*16)-1;
    ids_ignore{2,i} = flip(idstart(i+1)-(1:3*16));
end

for i = 1:nvert
    vstruc = chunkobj.vstruc{i};
    for j = 1:length(vstruc{1})
        jedge = vstruc{1}(j);
        for k = 1:length(vstruc{1})
        kedge = vstruc{1}(k);
        sysmat0(ids_ignore{1.5+vstruc{2}(j)/2,jedge}, ids_rcip{1.5+vstruc{2}(k)/2,kedge}) = 0;
        end
    end
end

val = sysmat0*dens;

ids = [];
fluxes = zeros(1,nedge);
for i = 1:nedge
    ids = edgeids(chunkobj,i);
    fluxes(i) = sum(val(ids).*chunkobj.wts(ids).');
end

nverts = size(chunkobj.verts,2);

ndepth = 40;
for ivert = 1:nverts
    starindtmp = rcipsav{ivert}.starind;
    
    solhat = dens(starindtmp);
    [solhatinterpcell,srcinfocell,wtscell] = chnk.rcip.rhohatInterp(solhat,rcipsav{ivert},ndepth);

    solinterp = [];
    % solhatinterp2 = [];
    chnkrs = cell(1,length(srcinfocell));
    for j = 1:length(srcinfocell)
        solhat = []; solhat.val = solhatinterpcell{j}(:).';

        chnkloc = chunkobj.echnks(chunkobj.vstruc{ivert}{1}(j));
        if chunkobj.vstruc{ivert}{2}(j) == -1
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

    verts = [[0;0]];
    edge2verts = [[NaN;NaN],[NaN;NaN]];
    cgrphloc = chunkgraph(verts,edge2verts,chnkrs);

    clear fkernloc
    if (size(kern) == 1)
        fkernloc = kern;
    else
        fkernloc(length(srcinfocell),length(srcinfocell)) = kernel();
        fkernloc(:,:) = kerns(chunkobj.vstruc{ivert}{1},chunkobj.vstruc{ivert}{1});
    end 

    opts = []; opts.rcip = false; %opts.adaptive_correction = true;
    sysmatloc = chunkermat(cgrphloc,fkernloc,opts);

    val_loc = sysmatloc * solinterp;

    % val_loc = val_loc + solinterp;

    for j = 1:length(srcinfocell)
        ids = edgeids(cgrphloc,j);
        flux_loc = sum(val_loc(ids).*cgrphloc.wts(ids).');

        fluxes(chunkobj.vstruc{ivert}{1}(j)) = ... 
            fluxes(chunkobj.vstruc{ivert}{1}(j)) + flux_loc;
    end
end

if isfield(opts,'tot')
    if opts.tot
        fluxes = sum(fluxes);
    end
end
end