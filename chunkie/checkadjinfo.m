function ier = checkadjinfo(chnkr)

i1 = 1;

hit = zeros(chnkr.nch,1);

for i = 1:chnkr.nch
    i2 = chnkr.adj(1,i1);
    hit(i2) = hit(i2)+1;
    i1 = i2;
end

ier = 0;
if nnz(hit == 0) > 0
    ier = ier + 1;
end
if nnz(hit > 1) > 0
    ier = ier + 2;
end
% fprintf('nch %d\n',chnkr.nch)
% fprintf('number hit %d\n',nnz(hit > 0))
% fprintf('number missed %d\n',nnz(hit == 0))
% fprintf('number doubled up %d\n',nnz(hit > 1))