%%%%%%%%%   KERNEL FOR LAYER POTENTIALS FOR OSCILLATORY BIHARMONIC %%%%%%


function submat= kern(k,srcinfo,targinfo,type,varargin)
  
src = srcinfo.r(:,:);
targ = targinfo.r(:,:);


% single layer 
if strcmpi(type,'s')
  submat = chnk.bihar2d.green(k,src,targ);
end

