function [WfrakL,accept,na]=WfrakLinit(A,trans,mscale,T)
% *** T is target vector, sources on canonical interval ***
ngl=length(T);
accept=1:ngl;
T=trans+mscale*T;
accept=accept(abs(T)<2);
na=length(accept);
p=zeros(ngl+1,1);
Q=zeros(ngl,na);
c=(1-(-1).^(1:ngl))./(1:ngl);
for j=1:na
  jj=accept(j);
  p(1)=log(abs((1-T(jj))/(1+T(jj))));
  p111=log(abs(1-T(jj)^2));
  for k=1:ngl
    p(k+1)=T(jj)*p(k)+c(k);
  end
  Q(1:2:ngl,j)=(p111-p(2:2:ngl))./(1:2:ngl)';
  Q(2:2:ngl,j)=(p(1)-p(3:2:ngl+1))./(2:2:ngl)';
end
WfrakL=Q.'/A;