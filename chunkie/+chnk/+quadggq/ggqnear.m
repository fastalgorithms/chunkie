function [xs1,ws1] = ggqnear(kmax)

switch kmax
    case 16
        [xs1,ws1] = chnk.quadggq.ggqnear16();
    case 20
        [xs1,ws1] = chnk.quadggq.ggqnear20();
    case 24
        [xs1,ws1] = chnk.quadggq.ggqnear24();
    case 30
        [xs1,ws1] = chnk.quadggq.ggqnear30();   
    case 40
        [xs1,ws1] = chnk.quadggq.ggqnear40();    
    case 60
        [xs1,ws1] = chnk.quadggq.ggqnear60();            
    otherwise
        error('near kmax selection not available')
end
