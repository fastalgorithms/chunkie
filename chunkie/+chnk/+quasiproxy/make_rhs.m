function rhs = make_rhs(chnkr,kh1,theta)

    [u_inc,un_inc] = make_incident(theta,kh1,chnkr);
    rhs=[-u_inc;-un_inc];

return

function [u_inc,un_inc] = make_incident(theta,kh,C)
% - Incident wave generator (for the top interface)
% INPUT: 
%       theta, incident angle
%       kh, wave number of the top layer
%       C, the parameterized geometry for the top interface

ima = sqrt(-1);
nn1 = C.n(1,:).';
nn2 = C.n(2,:).';


 u_inc = exp(ima*kh*(cos(theta)*C.r(1,:)'+sin(theta)*C.r(2,:)'));
 
 un_inc = ima*kh*(nn1*cos(theta)+nn2*sin(theta)).*u_inc;


return

