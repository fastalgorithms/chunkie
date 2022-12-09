  function [c,s]=rotmatc(a,b)
  if  a==0
    c=0;
    s=1;
  else
    temp=b/a;
    c=1/sqrt(1+abs(temp)^2);
    s=temp*c;
  end