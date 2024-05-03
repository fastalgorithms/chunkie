function intcoeffs = intpol(coeffs,const_option)
%LEGE.INTPOL compute the coefficients of the indefinite 
% integral of the polynomial with the coefficients given on 
% input
%
% input:
%  coeffs - array of Legendre series coefficients. if matrix, each column
%    contains the Legendres series coefficients of a separate function
% 
% optional input:
%  const_option - string. default ('true'). If 'true', then the constant on
%    output is such that the polynomial(s) described by intcoeffs on output
%    take the value zero at -1. If 'original', then the constant on output
%    is such that the polynomial(s) described by intcoeffs(1:end-1,:) on
%    output take the value zero at -1. This is used for a spectral
%    differentiation matrix taking point values to point values on a grid
%    of the same order.
%
% output:
%  intcoeffs - array of Legendre series coefficients for a definte integral
%    of the input Legendre series 
%

if nargin < 2
    const_option = 'true';
end

sz = size(coeffs);
szint = sz; szint(1) = szint(1)+1;
intcoeffs = zeros(szint);

if strcmpi(const_option,'true')
    ncc = szint(1);
elseif strcmpi(const_option,'original')
    ncc = szint(1)-1;
else
    error('LEGE.INTPOL: unknown option for constant');
end

for i = 2:sz(1)
    j = i-1;
    intcoeffs(i+1,:) = coeffs(i,:)/(2*j+1);
    intcoeffs(i-1,:) = -coeffs(i,:)/(2*j+1)+intcoeffs(i-1,:);
end

intcoeffs(2,:) = coeffs(1,:)+intcoeffs(2,:);

sss=-1;
dd = zeros(size(intcoeffs(end,:)));
for k = 2:ncc
    dd=dd+intcoeffs(k,:)*sss;
    sss=-sss;
end
intcoeffs(1,:)=-dd;

end