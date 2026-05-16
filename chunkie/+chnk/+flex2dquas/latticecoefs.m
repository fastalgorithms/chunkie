function S = latticecoefs(n,zk,d,kappa,alpha,a,M,l)
%CHNK.FLAX2DQUAS.LATTICECOEFS precompute lattice sum integrals 
% see CHNK.HELM2DQUAS.LATTICECOEFS

sn1 = chnk.helm2dquas.latticecoefs(n,zk,d,kappa,alpha,a,M,l);
sn2 = chnk.helm2dquas.latticecoefs(n,1i*zk,d,kappa,alpha,a,M,l);
S = cat(3,sn1,sn2);

end