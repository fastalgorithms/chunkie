function tau = taus(trap)

d = trap.d;

dd = sqrt(sum(abs(d).^2,1));
tau = bsxfun(@rdivide,trap.d,dd);

end

