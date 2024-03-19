% parameters

L = 20;

D = 12;

gam = 1.4;

bet = 1.2;

dx = 0.5;

Vr = 0.3;

len = 0:dx:L;

Kx = 300;

 

% initial conditions

P0 = 1;

V0 = D^2*(pi/4)*dx;

T0 = 25;

roh0 = 1.293;

Nu0 = 25;

R = 287;

tan_theta = 0.15;

 

% IGV calc

Eff = 100;