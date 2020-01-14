%SIMULATION INITIALIZATION Junior_Project_Sim 

FrontPlateLength = .130; %m
BackPlateLength = .260; %m
Width = .110; %m
rho = 1.225; %kg/m^3
v = 53.65; %m/s

q = 1/2*(rho)*(v)^2;

a = [0:0.0001:pi];
x = FrontPlateLength + BackPlateLength - ( BackPlateLength.*sqrt( 1- (FrontPlateLength/BackPlateLength)^2.*(sin(a)).^2 ) + FrontPlateLength.*cos(a) );

Ax = [x',a'];