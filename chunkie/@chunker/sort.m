function [obj,ifclosed] = sort(obj)
%SORT sort the chunker object so that adjacent chunks have sequential 
% indices
% 

ileft = find(obj.adj(1,:) < 1);
iright = find(obj.adj(2,:) < 1);

assert(length(ileft) <= 1,...
    'found more than one left end');
assert(length(iright) <= 1,...
    'found more than one right end');

ifclosed=true;
icurrent = 1;

if or(length(ileft) == 1,length(iright)==1)
    assert(and(length(iright)==1,length(ileft)==1),'chunker has one free end');
    icurrent = ileft;
    ifclosed=false;
end

adj2 = obj.adj(2,:);
inds = 1:obj.nch;
for i = 1:obj.nch
    inds(i) = icurrent;
    icurrent = adj2(icurrent);
end

assert(all(inds >= 1),'problem with found chunk indices');
if (length(unique(inds)) ~= length(inds))
    warning('sort behavior ill defined for chunkers with disjoint parts');
end

% reorder

obj.r = obj.r(:,:,inds);
obj.d = obj.d(:,:,inds);
obj.d2 = obj.d2(:,:,inds);
obj.h = obj.h(inds);

if obj.hasdata
    obj.data = obj.data(:,:,inds);
end

inds = 1:obj.nch;
obj.adj(1,:) = inds-1;
obj.adj(2,:) = inds+1;

obj.adj(1,1) = obj.nch;
obj.adj(2,obj.nch) = 1;

if ~ifclosed
    obj.adj(1,1) = -1;
    obj.adj(2,obj.nch) = -1;
end

    
end