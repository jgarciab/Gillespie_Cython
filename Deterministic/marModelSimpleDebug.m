

% dy/dt = b -  a*y
% dz/dt = g - a*z

function [dMARdt]=marModelSimpleDebug(t,MAR,la, lr, k_a, ka, k_r, kr, a00, a01, a10, lmar, br, ba)
P00 = MAR(1);
P01 = MAR(2);
marRA = MAR(3);
R2 = MAR(4);

dMARdt(1,1) = P01*k_r - P00*kr*R2;
dMARdt(2,1) = P00*kr*R2 - P01*k_r;
dMARdt(3,1) = P00*a00 + P01*a01 - marRA*lmar;
dMARdt(4,1) = marRA*br + P01*k_r  - P00*kr*R2  - R2*lr;
