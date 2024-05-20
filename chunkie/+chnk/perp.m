function n = perp(tau)

n = [tau(2,:); -tau(1,:)];
n = reshape(n,size(tau));