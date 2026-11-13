function [q0,q1,q0d] = qleg_half(t)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   caveat utilitor: this function evaluates Q_{-1/2} (and its
%    derivative) and Q_{1/2} at t+1;
%
% Input:
%     t - value between 0 \to \infty, corresponds to actual argument
%        of t+1 to Q funs
%
% Output:
%     q0 - value of Q_{-1/2}
%     q1 - value of Q_{1/2}
%     q0d - value of Q'_{-1/2}
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    inear = find(t<0.01);
    ifar  = find(t>100);
    imid  = find((t>=0.01).*(t<=100));
    [qmm,qmpm] = chnk.axissymlap2d.qget_zero(1+t(imid));
    [qmn] = chnk.axissymlap2d.qeval0_near(t(inear));
    [qmpn] = chnk.axissymlap2d.qeval1_near(t(inear));
    [qmdn] = chnk.axissymlap2d.qeval0der_near(t(inear));
    [qmf] = chnk.axissymlap2d.qeval0_far(1+t(ifar));
    [qmpf] = chnk.axissymlap2d.qeval1_far(1+t(ifar));


    
    q0 = zeros(size(t));
    q0(imid)  = qmm;
    q0(inear) = qmn;
    q0(ifar)  = qmf;

    q1 = zeros(size(t));
    q1(imid)  = qmpm;
    q1(inear) = qmpn;
    q1(ifar)  = qmpf;

    q0d = -(1/2)*q1 + (1/2)*(1+t).*q0;
    q0d = q0d./(-t.*(t+2));
    q0d(inear) = qmdn;
end
