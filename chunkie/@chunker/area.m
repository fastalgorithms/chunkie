function a = area(chnkr)

wts = whts(chnkr);
rnorm = normals(chnkr);
a = sum(sum(wts.*squeeze(sum(rnorm.*chnkr.r,1))))/3.0;

end