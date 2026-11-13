function [qm,qmp] = qget_zero(x)
%
% chnk.axissymlap2d.qget_zero evaluates the 
%
        a = sqrt(2./(x+1));
        [fF,fE] = chnk.axissymlap2d.gaus_agm(x);

        qm = a.*fF;
        qmp = x.*qm-2*fE./(a);

end
