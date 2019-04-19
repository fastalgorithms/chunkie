function rmax = max(obj)

rmax = max(reshape(obj.r,obj.dim,obj.k*obj.nch),[],2);