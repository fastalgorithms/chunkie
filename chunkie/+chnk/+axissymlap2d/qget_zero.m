function [qm,qmp] = qget_zero(x)

        a = sqrt(2./(x+1));
        [fF,fE] = chnk.axissymlap2d.gaus_agm(x);

        q0 = a.*fF;
        q1 = x.*q0-2*fE./(a);

        qa = q0;
        qb = q1;

        qmm = 0;
        qm  = q0;
        qmp = q1;

end