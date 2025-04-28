function chnkr = chunkerpoints(r)

pref = [];

ifclosed=true;
dim = size(r,1);
k   = size(r,2);
nch = size(r,3); 
 
pref.dim = dim;
pref.k   = k;

[xs,~,us,vs] = lege.exps(k);
dermat = (vs*[lege.derpol(us); zeros(1,k)]).';


%       . . . start chunking

adjs = zeros(2,nch);

adjs(1,:) = (1:nch) - 1;
adjs(2,:) = (1:nch) + 1;

adjs(1,1)=nch;
adjs(2,nch)=1;


%       up to here, everything has been done in parameter space, [ta,tb]
%       . . . finally evaluate the k nodes on each chunk, along with 
%       derivatives and chunk lengths

chnkr = chunker(pref); % empty chunker
chnkr = chnkr.addchunk(nch);


for i = 1:nch
    
    rtmp  = squeeze(r(:,:,i));
    dtmp  = rtmp*dermat;
    d2tmp = dtmp*dermat;
    chnkr.rstor(:,:,i) = rtmp;
    chnkr.dstor(:,:,i) = dtmp;
    chnkr.d2stor(:,:,i) = d2tmp;
end

chnkr.adjstor(:,1:nch) = adjs(:,1:nch);

% Set normals
chnkr.nstor(:,:,1:nch) = normals(chnkr);

% Set weights
chnkr.wtsstor(:,1:nch) = weights(chnkr);

end
