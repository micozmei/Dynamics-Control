%% INITIALIZE VAWT SIM
clc; close all; clear all;

%% Misc Parameters
rho = 1.225; %fluid (air) density; [kg/m^3]

%% Airfoil Numbers
% Values estimated from
A0 = 0;             % Alpha where C_l = 0; [degree]
CL1_max = 1.4;      % max C_L pre-stall; [dimentionless]
CL2_max = 1.2;      % max C_L post-stall; [dimentionless]
S1 = 0.6/6;         % slope of C_l vs alpha graph at alpha = A0; [1/degrees]
CD0 = 0.005;        % CD at A0; [dimentionless]
CD1_max = 0.017;    % CD at ACD1; [dimentionless]
CD2_max = 1.8;     % CD at ACD2 (90 degrees); [dimentionless]
ACD1 = 12;          % alpha at stall transition; [degree]
ACL1 = 15;          % alpha at stall; [degree]

%ACD1up = 360 - ACD1;%[degree]


%% CL NACA0012 Lookup Table values
angles =       [ 0   10   15   20   22    30     45     60    75   90   104    120    135   150   160    165   172    180  188    195  200    210  225  240 256   270 285   300     315 330 338 340 345 350 360];
C_l_NACA0012 = [ 0   1.15 1.4  0.9  0.91  1.05   1.1    0.95  0.6  0.0  -0.44  -0.82  -1.1  -0.9  -0.75  -0.8  -0.82  0   0.82   0.8   0.75   0.9  1.1  0.8  0.44  0.0  -0.6  -0.95  -1.1  -1.05 -0.91 -0.9  -1.4  -1.15  0 ];
C_l_NACA0012_spline = spline(angles, C_l_NACA0012, [0:1:360]);

% plot to compare with reference graph
% figure(2)
% plot(angles, C_l_NACA0012,'o',[0:1:360],C_l_NACA0012_spline)

%% Blade Parameters

Blade_Length = 40;                   %turbine blade length; [m] 
Chord = 4;                      %turbine blade chord; [m]
Planform_A = Blade_Length*Chord;    %planfrom blade area; [m^2]
AR = Blade_Length/Chord;            %aspect ratio; [dimentionless]
e = 0.70;                           %oswald efficiency factor, rectangular wing

%% Generator Parameters
% c_wg = 0.000290; %generator linear friction coefficient
% c_wr = 0.000290; %rotor linear friction coefficient
% %Cg = .001; %Generator damping
% %C_Tg = .001; %Generator torsional damping
% %Cr = .001; %Rotor damping 
% %C_Tr = .001;%Rotor torsional damping
% eta_g = 0.80; %Generator efficiency
% J_g = ((0.01)^4)*(pi/2); %Generator rotational inertia
% K = 0.084; %Generator machine constant
% Kr = ((27*10^9)*((pi*(0.01)^4)/2))/(.1); %Rotor torsional stiffness
% Kg = ((27*10^9)*((pi*(0.01)^4)/2))/(.1); %Generator torsional stiffness
% L = 0.0023; %Generator inductance
% N = 2.3; %Gear Ratio
% Re = 0.05; %Monitoring resistor?
% R2 = 1.51; %Generator armature resistance
% tau_e = 1;
% Bg = 0.001;

Jw = 1.6*10^6;
Jg = 15.8;
K = 6*10^7;
D = 10^6;
N = 100;
eta_g = 0.9;

%Eigenswings
a11 = 0.8;
a12 = 0.2;
a21 = 0.5;
a22 = 0.5;
phi11 = 0;
phi12 = pi/2;
phi21 = 0;
phi22 = pi/2;
h1 = 1;
h2 = 1;
a1 = 0.01;
a2 = 0.08;

%wind variation
A_wind = 10;
w_wind = (2*pi)*5;
v_wind = 5;
w_gust = (2*pi)/(30);
v_gustmax = v_wind;