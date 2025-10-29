function [cgrph, isort, rfins, matkerns, plotkerns, idgnore] = build_fins(cgrph, matkerns, plotkerns, iverts, iedges, iedges_ref, nchs, cs)
% Wrapper to add fins to a chunkgraph and update the kernels
% 
% maintains number of discretization nodes
%
% Input: 
%   cgrph - chunkgraph
%   matkerns - matrix of kernels for the chunkermat 
%   plotkerns - vector of kernels for chunkerkerneval
%   iverts - vertex numbers
%   iedges - edge of verts to rotate 
%   iedges_ref - edge to cancel
%   nchs - number of chunks to use
%   cs - rotate kernel coefficient
% Output:
%   cgrph - split chunkgraph
%   isort - mapping from initial nodes to new nodes
%   rfins - cell array of fins
%   idingore - list of extra fake vertices for rcip to ignore


% verify shape of kernels
nedge = length(cgrph.echnks);
assert(all(size(matkerns) == [nedge, nedge]))
assert(size(plotkerns,2) == nedge)

% number of fins
nfin = length(iverts);

% cell array of fins
rfins = cell(1,nfin);

% initial number of verts
nvert0 = size(cgrph.verts,2);

isort = 1:cgrph.npt;
% add each fin
for i = 1:nfin
    % get fin info
    ivert = iverts(i);
    iloc_edge = iedges(i);
    iedge = cgrph.vstruc{ivert}{1}(iloc_edge);

    nch_loc = cgrph.echnks(iedge).nch;
    if cgrph.vstruc{ivert}{2}(iloc_edge) == -1
        chsplit = nchs(i);
        irefl = iedge;
        inrefl = length(cgrph.echnks)+1;

        [~,tau_s] = chunkends(cgrph.echnks(iedge), 1);
        tau_s = tau_s(:,1);
    else
        chsplit = nch_loc - nchs(i);
        inrefl = iedge;
        irefl = length(cgrph.echnks)+1;

        [~,tau_s] = chunkends(cgrph.echnks(iedge), nch_loc);
        tau_s = tau_s(:,2);
    end
    iedge_t =  cgrph.vstruc{ivert}{1}(iedges_ref(i));
    if cgrph.vstruc{ivert}{2}(iedges_ref(i)) == -1
        [~,tau_t] = chunkends(cgrph.echnks(iedge_t), 1);
        tau_t = tau_t(:,1);
    else
        [~,tau_t] = chunkends(cgrph.echnks(iedge_t), cgrph.echnks(iedge_t).nch);
        tau_t = tau_t(:,2);
    end

    % split chunkgraph
    [cgrph,isort_tmp] = chnk.fins.split_chunkgraph(cgrph, iedge, chsplit);
    isort = isort(isort_tmp);

    % update kernels
    nedge_new = length(cgrph.echnks);
    matkerns(nedge_new,:) = matkerns(iedge,:);
    matkerns(:,inrefl) = matkerns(:,iedge);
    plotkerns(:,inrefl) = plotkerns(:,iedge);

    % rotation parameters
    r_vert = cgrph.verts(:,ivert);
    alpha = -2 * (atan2(tau_s(2),tau_s(1)) - atan2(tau_t(2),tau_t(1)) );   
    for j = 1:nedge_new
        % TODO: add reflections
        matkerns(j,irefl) = matkerns(j,inrefl);
        plotkerns(:,irefl) = plotkerns(:,inrefl);
    end

    % make new fin
    rot = [cos(alpha), -sin(alpha); sin(alpha), cos(alpha)];
    rfins{i} = chnk.fins.rotate_pt(cgrph.echnks(irefl),r_vert, rot);
end

idgnore = (nvert0+1):size(cgrph.verts,2);
end

