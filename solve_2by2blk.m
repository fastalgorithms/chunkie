function [sol] = solve_2by2blk(rhs,Fs1,Fs2,skel,opts_perm)
    iperm = opts_perm.iperm;
    rhs   = rhs(iperm,:);
    rhsu  = rhs(1:opts_perm.ns(1),:);
    rhsl  = rhs((opts_perm.ns(1)+1):end,:);
    
    sol11 = rskelf_sv(Fs1,rhsu)+skel.U11*(skel.V11*rhsu);
    sol12 = skel.U12*(skel.V12*rhsl);
    sol21 = skel.U21*(skel.V21*rhsu);
    sol22 = rskelf_sv(Fs2,rhsl)+skel.U22*(skel.V22*rhsl);
    sol = [sol11+sol12;sol21+sol22];
    sol = sol(opts_perm.invperm);
end

