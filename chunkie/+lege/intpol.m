function intcoeffs = intpol(coeffs)
%LEGE.INTPOL compute the coefficients of the indefinite 
% integral of the polynomial with the coefficients given on 
% input

sz = size(coeffs);
szint = sz; szint(1) = szint(1)+1;
intcoeffs = zeros(szint);

for i = 2:sz(1)
    j = i-1;
    intcoeffs(i+1,:) = coeffs(i,:)/(2*j+1);
    intcoeffs(i-1,:) = -coeffs(i,:)/(2*j+1)+intcoeffs(i-1,:);
end

intcoeffs(2,:) = coeffs(1,:)+intcoeffs(2,:);

sss=-1;
dd = zeros(size(intcoeffs(end,:)));
for k = 2:(szint(1)-1)
    dd=dd+intcoeffs(k,:)*sss;
    sss=-sss;
end
intcoeffs(1,:)=-dd;

end