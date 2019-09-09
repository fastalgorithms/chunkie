function rmin = min(obj)

rmin = min(reshape(obj.r,obj.dim,obj.k*obj.nch),[],2);