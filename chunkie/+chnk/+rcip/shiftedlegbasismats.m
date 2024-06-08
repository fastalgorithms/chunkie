function [sbclmat,sbcrmat,leftvalmat,rightvalmat,u] = shiftedlegbasismats(k)
%shiftedlegbasismats
%
% utility for rcip. not designed for end-user 
% 
% input 
% k - legendre node order on chunk 
%
% output
% sbclmat - matrix from function values at legendre
%                nodes scaled to [0,1] to expansion 
%                of the form t*\sum_{j=0}^{k-1} c_j P_j(2t-1) 
% sbcrmat - matrix from function values at legendre
%                nodes scaled to [-1,0] to expansion 
%                of the form t*\sum_{j=0}^{k-1} c_j P_j(2t+1) 
% leftvalmat - matrix from function values to value of 
%              interpolant at left end
%

% author: Travis Askham
  
[t0,~,u,v] = lege.exps(k);

t = (t0+1)/2;
basis = t(:).*v(:,1:end-1);
[uu,ss,vv] = svd(basis);

sbclmat = vv(:,1:k-1)*(ss(1:k-1,1:k-1)\(uu(:,1:k-1)'));

t = (t0-1)/2;
basis = t(:).*v(:,1:end-1);
[uu,ss,vv] = svd(basis);

sbcrmat = vv(:,1:k-1)*(ss(1:k-1,1:k-1)\(uu(:,1:k-1)'));

pm1 = lege.pols(-1,k-1);
leftvalmat = (pm1(:).')*u;

p1 = lege.pols(1,k-1);
rightvalmat = (p1(:).')*u;

end
