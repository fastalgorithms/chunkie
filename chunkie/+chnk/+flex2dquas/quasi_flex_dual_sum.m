function [val, grad, hess, third, fourth] = quasi_flex_dual_sum(rx,ry,zk,kappa,d)

[valh0, gradh0, hessh0, thirdh0, fourthh0] = chnk.lap2dquas.quasi_dual_sum(rx,ry,zk,kappa,d);
[valk0, gradk0, hessk0, thirdk0, fourthk0] = chnk.lap2dquas.quasi_dual_sum(rx,ry,1i*zk,kappa,d);

over2k2 = 1; %/(2*zk.^2);
val = over2k2*(valh0-valk0);
grad = over2k2*(gradh0-gradk0);
hess = over2k2*(hessh0-hessk0);
third = over2k2*(thirdh0 - thirdk0);
fourth = over2k2*(fourthh0-fourthk0);

end