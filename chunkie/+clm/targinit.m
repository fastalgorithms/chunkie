function [xg,yg,targs,ngrtot]=targinit(xylim,ngr)
xg=linspace(xylim(1),xylim(2),ngr);
yg=linspace(xylim(3),xylim(4),ngr);
ngrtot=ngr^2;
targs=zeros(2,ngrtot);
for k=1:ngr
  targs(1,(k-1)*ngr+(1:ngr)) = xg(k);
  targs(2,(k-1)*ngr+(1:ngr)) = yg; 
end
