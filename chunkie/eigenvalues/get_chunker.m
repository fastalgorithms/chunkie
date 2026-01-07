function [chnkr, nholes] = get_chunker(j)
rad = 1; ctr = [0.0;0.0];
circfun = @(t) ctr + rad*[cos(t(:).');sin(t(:).')];

opts = [];
opts.maxchunklen = 0.5;
pref = []; pref.nchmax = 100000;
chnkr1 = chunkerfunc(circfun, opts, pref);
plot(chnkr1, 'k.');

l = j;
chnkreps = chnkr1;
nh = 0;
k = 1/(4 * l); 
h = 5 * l;
chunkerhole0 = k * chnkr1;
chunkerhole0 = chunkerhole0.reverse();
for i = 0 : l-1
    for j = 0 : l-1
        check1 = [i/l,j/l];
        l1 = norm(check1 + [1/h,1/h]); 
        l2 = norm(check1 + [1/l,0.0] + [-1/h,1/h]);
        l3 = norm(check1 + [0.0,1/l] + [1/h,-1/h]);
        l4 = norm(check1 + [1/l,1/l] + [-1/h,-1/h]);
        if (l1<=1) && (l2<=1) && (l3<=1) && (l4<=1)
            chnkrctr1 = check1 + (1/(2*l))*[1.0,1.0];
            chnkrctr2 = [-1,1].*chnkrctr1;
            chnkrctr3 = [1,-1].*chnkrctr1;
            chnkrctr4 = [-1,-1].*chnkrctr1;
            chnkrhole1 = chnkrctr1 + chunkerhole0;
            chnkrhole2 = chnkrctr2 + chunkerhole0;
            chnkrhole3 = chnkrctr3 + chunkerhole0;
            chnkrhole4 = chnkrctr4 + chunkerhole0;
            chnkreps = merge([chnkreps, chnkrhole1, chnkrhole2, chnkrhole3, chnkrhole4]);
            nh = nh + 4;
         end
    end
 end
    chnkr = chnkreps;
    nholes = nh;
end