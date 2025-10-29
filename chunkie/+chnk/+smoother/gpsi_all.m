function [val,grad,hess,hess_sig] = gpsi_all(dx,dy,dr,dsigt)

    [val,vald] = chnk.smoother.psi_eval(dr,dsigt);

    da = 1/(4*pi);
    z  = dr.^2./(2*dsigt.*dsigt);
    val = val + log(z)*da;

    if (nargout==1) 
        return
    end

    grad = zeros(numel(dx),2);
    grad(:,1) = vald.*dx./dr;
    grad(:,2) = vald.*dy./dr;

    grad(:,1) = grad(:,1)+2*dx./(dr.*dr)*da;
    grad(:,2) = grad(:,2)+2*dy./(dr.*dr)*da;

    if (nargout==2)
        return
    end

    hess = zeros(numel(dx),2,2);

    hess(:,1,1) = 0;
    hess(:,1,2) = 0;
    hess(:,2,1) = 0;
    hess(:,2,2) = 0;

    eterm = exp(-dr.*dr./(2*dsigt.^2));

    hess(:,1,1) = eterm./(dr.*dr);
    hess(:,2,2) = eterm./(dr.*dr);

    fact = -(2*dsigt.^2+dr.^2)./(dsigt.^2.*dr.^2).*eterm;
    hess(:,1,1) = hess(:,1,1) + (dx.*dx./dr.^2).*fact;
    hess(:,1,2) = hess(:,1,2) + (dx.*dy./dr.^2).*fact;
    hess(:,2,1) = hess(:,2,1) + (dx.*dy./dr.^2).*fact;
    hess(:,2,2) = hess(:,2,2) + (dy.*dy./dr.^2).*fact;

    hess(:,1,1) = hess(:,1,1)/(2*pi);
    hess(:,2,1) = hess(:,2,1)/(2*pi);
    hess(:,1,2) = hess(:,1,2)/(2*pi);
    hess(:,2,2) = hess(:,2,2)/(2*pi);

    hess(:,1,1) = hess(:,1,1) - 2./(dr.*dr)*da;
    hess(:,2,2) = hess(:,2,2) - 2./(dr.*dr)*da;

    fact = da*4./dr.^4;
    hess(:,1,1) = hess(:,1,1) + dx.*dx.*fact;
    hess(:,1,2) = hess(:,1,2) + dx.*dy.*fact;
    hess(:,2,1) = hess(:,2,1) + dx.*dy.*fact;
    hess(:,2,2) = hess(:,2,2) + dy.*dy.*fact;

    if (nargout==3)
        return
    end


    fact = eterm./(dsigt.^3)/(2*pi);
    hess_sig(:,1) = -dx.*fact;
    hess_sig(:,2) = -dy.*fact;

end
