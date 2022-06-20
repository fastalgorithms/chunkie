function [nearinfo,varargout] = interp2rel(chnkinfo,intstruc)

    if (~isfield(chnkinfo,'r'))
        error('In interp2rel, chnkinfo must have an "r" field');
    end
    if (~isfield(chnkinfo,'d'))
        error('In interp2rel, chnkinfo must have an "d" field');
    end
    if (~isfield(chnkinfo,'d2'))
        error('In interp2rel, chnkinfo must have an "d2" field');
    end
    
    n = size(chnkinfo.r,2);
    
    if (nargin <2 || ~isfield(intstruc,'v2cmat') || ... 
       ~isfield(intstruc,'vcint') )
        
        [xl,wl,ul,vl]     = lege.exps(n);
        [xl2,wl2,ul2,vl2] = lege.exps(n-1);

        [pols,~] = lege.pols(-1,n-1);
        [vfns,~] = lege.pols(xl2,n-1);

        vlmod  = vfns' - ones([n-1,1])*pols';
        wei_mat = 2./(1+xl2)*ones([1,n]);
        vlmod(:,1) = 0;

        v2cmat = ul2*(vlmod.*wei_mat)*ul;
        
        intstruc = [];
        intstruc.v2cmat = v2cmat;
        
        [aint,~,~] = lege.intmat(n);
        intstruc.vcint = v2cmat*aint;
        
        [pols,ders] = lege.pols(xl2,n-1);
        pols = pols'*ul;
        ders = ders'*ul;
        intstruc.vvmat = pols;
        intstruc.vdmat = ders;
        intstruc.vl    = vl2;
        ps = lege.pols(xl,n-2);
        intstruc.eval  = ps';
    end
    
    nearinfo = [];
    
    dx = chnkinfo.d(1,:);
	dy = chnkinfo.d(2,:);
    
    x   = intstruc.vcint*dx';
	y   = intstruc.vcint*dy';
	nearinfo.r = [x';y'];
    
    xx  = chnkinfo.r(1,:);
    yy  = chnkinfo.r(2,:);
    %x   = intstruc.v2cmat*xx';
    %y   = intstruc.v2cmat*yy';
    %nearinfo.r = [x';y'];
    
	dxt = intstruc.vvmat*dx';
	dyt = intstruc.vvmat*dy';

    d = [dxt',dyt'];
    
    dx2 = chnkinfo.d2(1,:)';
	dy2 = chnkinfo.d2(2,:)';
    
    dxt2 = intstruc.vvmat*dx2;
	dyt2 = intstruc.vvmat*dy2;
    
    d2 = [dxt2';dyt2'];

    nearinfo.d  = d;
    nearinfo.d2 = d2;
    
    nearinfo.d  = chnkinfo.d;
    nearinfo.d2 = chnkinfo.d2;
    
    if (isfield(chnkinfo,'data'))
        nearinfo.data = chnkinfo.data;
    end
    
    if (nargout == 2)
        varargout{1} = intstruc;
    end
    
end

