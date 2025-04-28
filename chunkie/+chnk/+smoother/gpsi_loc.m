function [val,grad,hess,hess_sig] = gpsi_loc(dx,dy,dr,dsigt)

    degam = 0.57721566490153286060651209008240243104;
    da = sqrt(2)*dsigt;
    z = dr./da;
    z2 = z.*z;
    val = (z2-z2.^2/4+z2.^3/18-z2.^4/96+z2.^5/600-z2.^6/4320);
    val = val./(4*pi);
    val = val - degam/(4*pi); 

    if (nargout==1) 
        return
    end

    grad = zeros(numel(dx),2);
    vvd = (1-z2/2+z2.^2/6-z2.^3/24+z2.^4/120-z2.^5/720);
    grad(:,1) = dx.*vvd./(2*pi*da.^2);
    grad(:,2) = dy.*vvd./(2*pi*da.^2);

    if (nargout==2) 
        return
    end

    hess = zeros(numel(dx),2,2);
    hess(:,1,1) = -vvd./(2*pi*da.^2);
    hess(:,2,2) = -vvd./(2*pi*da.^2);
    hess(:,1,2) = 0;
    hess(:,2,1) = 0;

    vvvd = (-1+4*z2/6-6*z2.^2/24+8*z2.^3/120-10*z2.^4/720);
    vvvd = -vvvd./(2*pi*da.^2)./da./da;
    hess(:,1,1) = hess(:,1,1) + dx.*dx.*vvvd;
    hess(:,2,2) = hess(:,2,2) + dy.*dy.*vvvd;
    hess(:,1,2) = hess(:,1,2) + dx.*dy.*vvvd;
    hess(:,2,1) = hess(:,2,1) + dx.*dy.*vvvd;

    if (nargout==3)
        return
    end

    hess_sig = zeros(numel(dx),2);

    eterm = exp(-dr.*dr./(2*dsigt.^2));
    fact = eterm./(dsigt.^3)/(2*pi);
    hess_sig(:,1) = -dx.*fact;
    hess_sig(:,2) = -dy.*fact;

end
