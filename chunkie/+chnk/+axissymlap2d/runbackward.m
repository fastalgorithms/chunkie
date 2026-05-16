function [q0,q1,qd0,qdd0] = runbackward(xm1,n)
    % returns q0 = Q_{n/2-2}
    if rem(n, 2) == 0
        % n is even
        first_term = 1/2*log((xm1+2)./xm1); % Q0
        m = n/2-1;
    else
        % n is odd, n>1
        first_term = chnk.axissymlap2d.qleg_half(xm1); % Q_{-1/2}
        m = (n-1)/2;
    end
    j = 0;
    jmax = 300;
    tol = 1e-8; 
    t = repmat(zeros(size(xm1)),1,1,2);
    t(:,:,1) = ones(size(xm1));
    t(:,:,2) = zeros(size(xm1));
    onezero = t;
    told = repmat(zeros(size(xm1)),1,1,2);
    told(:,:,1) = realmax*ones(size(xm1));
    told(:,:,2) = ones(size(xm1));
    x = xm1+1;
    while j<jmax
        v = n/2-1+j; %m+j-1/2;
        qv1 = t(:,:,2);
        qv0 = t(:,:,1); 
        qv2 = (2*v+3)/(v+2)*x.*qv1-(v+1)/(v+2)*qv0;
        t(:,:,1) = qv1;
        t(:,:,2) = qv2;
        vecn = vecnorm(t,2,3);
        t(:,:,1) = t(:,:,1)./vecn;
        t(:,:,2) = t(:,:,2)./vecn;
        t = t./t(:,:,2); 
        if max(abs(told(:,:,1)-t(:,:,1)),[],'all')<tol
            break;
        end 
        told = t;
        j = j+1;
    end
    t = onezero;
    for ii = m+j:-1:m
        qv1 = t(:,:,1);
        qv2 = t(:,:,2);
        qv0 = (2*v+3)/(v+1)*x.*qv1-(v+2)/(v+1)*qv2;
        t(:,:,1) = qv0;
        t(:,:,2) = qv1;
        t = t./qv0;
        v = v-1;
    end
    qm1 = t(:,:,2);
    for ii = m-1:-1:1
        qv1 = t(:,:,1);
        qv2 = t(:,:,2);
        qv0 = (2*v+3)/(v+1)*x.*qv1-(v+2)/(v+1)*qv2;
        t(:,:,1) = qv0;
        t(:,:,2) = qv1;
        v = v-1;
    end
    
    scale = first_term./t(:,:,1);
    q0 = scale; 
    q1 = qm1.*scale; 
    v = n/2-2; %
    qd0 = (-v-1)*q1+(v+1)*(1+xm1).*q0; 
    qd0 = -qd0./(2+xm1)./xm1;
    qd1 = (v+1)*q0-(v+1)*(1+xm1).*q1;
    qd1 = -qd1./(2+xm1)./xm1;
    qdd0 = (-v-1)*qd1+(v+1)*q0+(v+3)*(1+xm1).*qd0;
    qdd0 = -qdd0./(2+xm1)./xm1;
end