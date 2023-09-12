function [chnkrup,sigmaup] = upsample(chnkr,kup,sigma)
%UPSAMPLE return a chunker object which has been upsampled to a finer grid
%
% [chnkrup,sigmaup] = upsample(chnkr,kup,sigma);
%
% input:
% chnkr - chunker object
% kup - integer, order of panels for upsampled chunker
% sigma - optional, function on chunker to upsample
%
% output:
%
% chnkrup - upsampled verson of chunker object 
% sigmaup - upsampled function (if sigma provided), empty otherwise
% 

k = chnkr.k;
dim = chnkr.dim;
nch = chnkr.nch;

[~,~,u] = lege.exps(k);
[tu,wu,~,vu] = lege.exps(kup);

upmat = vu(:,1:k)*u;

pref = []; pref.k = kup;
chnkrup = chunker(pref,tu,wu);
chnkrup = chnkrup.addchunk(nch);

nn = dim*nch;

chnkrup.r = permute(reshape(upmat*reshape( ...
    permute(chnkr.r,[2,1,3]),k,nn),kup,dim,nch),[2,1,3]);
chnkrup.d= permute(reshape(upmat*reshape( ...
    permute(chnkr.d,[2,1,3]),k,nn),kup,dim,nch),[2,1,3]);
chnkrup.d2= permute(reshape(upmat*reshape( ...
    permute(chnkr.d2,[2,1,3]),k,nn),kup,dim,nch),[2,1,3]);
chnkrup.adj= chnkr.adj;
chnkrup.h= chnkr.h;

chnkrup.n = normals(chnkrup);

if chnkr.hasdata
    ndata = chnkr.size(data,1);
    nndata = ndata*nch;
    chnkrup.r = permute(reshape(upmat*reshape( ...
        permute(chnkr.datastor,[2,1,3]),k,nndata),kup,ndata,nch),[2,1,3]);
end

if nargin > 2
    dimsig = numel(sigma)/(nch*k);
    sigma= reshape(sigma,dimsig,k,nch);
    nnsig = dimsig*nch;
    sigmaup = permute(reshape(upmat*reshape( ...
        permute(sigma,[2,1,3]),k,nnsig),kup,dimsig,nch),[2,1,3]);
else
    sigmaup = [];
end
