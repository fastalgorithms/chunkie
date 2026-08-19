function [q0,q1,qd0,qdd0] = runforward(xm1,n)
    % returns q0 = Q_{n/2-2}
    if rem(n, 2) == 0
        % n is even
        Q1 = 1/2*log((xm1+2)./xm1); % Q0
        Q2 = (xm1+1)/2.*log((xm1+2)./xm1)-1; % Q1
        m = n/2-1;
        v = 0;
    else
        % n is odd, n>1
        [Q1,Q2,~,~] = chnk.axissymlap2d.qleg_half(xm1); % Q_{-1/2}
        m = (n-1)/2;
        v = -1/2;
    end
    t(:,:,1) = Q1;
    t(:,:,2) = Q2;
    while v < m-3/2
        qv0 = t(:,:,1);
        qv1 = t(:,:,2);
        qv2 = (2*v+3)/(v+2)*(xm1+1).*qv1-(v+1)/(v+2)*qv0;
        t(:,:,1) = qv1; 
        t(:,:,2) = qv2;
        v = v+1;
    end
    q0 = t(:,:,1); % Q_{m-3/2}
    q1 = t(:,:,2); % Q_{m-1/2}
    qd0 = (-v-1)*q1+(v+1)*(1+xm1).*q0; 
    qd0 = -qd0./(2+xm1)./xm1;
    qd1 = (v+1)*q0-(v+1)*(1+xm1).*q1;
    qd1 = -qd1./(2+xm1)./xm1;
    qdd0 = (-v-1)*qd1+(v+1)*q0+(v+3)*(1+xm1).*qd0;
    qdd0 = -qdd0./(2+xm1)./xm1;
end
