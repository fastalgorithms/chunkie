function [cgrph, isort, rfins, idgnore] = build_fins(cgrph, kerns, iverts, iedges, iedges_ref, nchs, cs)
% fin_info an array of structs, 
%   * iverts = vertex numbers
%   * iedges = edge of verts to rotate 
%   * nchs = number of chunks to use
%   * cs = rotate kernel coefficient


nfin = length(iverts);
kerns_0 = kerns;

rfins = cell(1,nfin);

nvert0 = size(cgrph.verts,2);

isort = 1:cgrph.npt;
for i = 1:nfin
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

    [cgrph,isort_tmp] = chnk.fins.split_chunkgraph(cgrph, iedge, chsplit);
    isort = isort(isort_tmp);

    nedge_new = length(cgrph.echnks);
    kerns(nedge_new,:) = kerns(iedge,:);
    kerns(:,inrefl) = kerns(:,iedge);

    r_vert = cgrph.verts(:,ivert);
    alpha = -2 * (atan2(tau_s(2),tau_s(1)) - atan2(tau_t(2),tau_t(1)) );   
    for j = 1:nedge_new
        kerns(j,irefl) = kerns(j,iedge);
    end

    rot = [cos(alpha), -sin(alpha); sin(alpha), cos(alpha)];
    rfins{i} = chnk.fins.rotate_pt(cgrph.echnks(irefl),r_vert, rot);
end

idgnore = (nvert0+1):size(cgrph.verts,2);



end

