function A = build_v2emat(obj)
%V2EMAT get a "vertex to edge" matrix for the chunkgraph object
% 
% input: obj - a chunkgraph object 
% 
% output: A - a nedges x nverts sparse matrix 
%  A(i,j) = 1 if edge i starts at vertex j
%  A(i,j) = -1 if edge i ends at vertex j
%  A(i,j) = 2 if edge i starts and ends at vertex j
%

nverts = size(obj.verts,2);
nedges = size(obj.edgesendverts,2);


ii = zeros(2*nedges,1);
jj = zeros(2*nedges,1);
vv = zeros(2*nedges,1);
nf = 0;

for i = 1:nedges
    j1 = obj.edgesendverts(1,i);
    j2 = obj.edgesendverts(2,i);
    if ~isnan(j1) && ~isnan(j2)
        if j1 == j2
            nf = nf + 1;
            ii(nf) = i;
            jj(nf) = j1;
            vv(nf) = 2;
        else
            nf = nf + 1;
            ii(nf) = i;
            jj(nf) = j1;
            vv(nf) = -1;
            nf = nf + 1;
            ii(nf) = i;
            jj(nf) = j2;
            vv(nf) = 1;
        end
    end
end

A = sparse(ii(1:nf),jj(1:nf),vv(1:nf),nedges,nverts);

