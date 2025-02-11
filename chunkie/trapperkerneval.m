function fints = trapperkerneval(trap,kern,dens,targobj,opts)
%TRAPPERINTKERN compute the convolution of the integral kernel with

wts = trap.wts;
srcinfo = []; srcinfo.r = trap.r; srcinfo.d = trap.d; srcinfo.d2 = trap.d2;
srcinfo.n = trap.n;

% Assign appropriate object to targinfo
targinfo = [];
if isa(targobj, "trapper")
    targinfo.r = targobj.r(:,:);
    targinfo.d = targobj.d(:,:);
    targinfo.d2 = targobj.d2(:,:);
    targinfo.n = targobj.n(:,:);
else
    targinfo.r = targobj;
end

mat = kern(srcinfo,targinfo);
fints = mat*diag(wts)*dens;

end
