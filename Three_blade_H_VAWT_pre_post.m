%% VAWT SIM
clc; close all; clear all;

%% ------------PRE-PROCESSING: ------------ %%

% Misc Parameters

rho = 1.225; %fluid (air) density; [kg/m^3]
blade1_radius = 20; %distance from center to foil 1; [m]
blade2_radius = 20; %distance from center to foil 2; [m]
blade3_radius = 20; %distance from center to foil 3; [m]

blade1_pitch = 0; %pitch of foil 1; [degrees]
blade2_pitch = 0; %pitch of foil 2; [degrees]
blade3_pitch = 0; %pitch of foil 3; [degrees]

wind_vector_i = [50,0,0]; %inertial wind vector; [m/s]

%angular_velocity_init = 1.048; %initial phi_dot; [rad/s]


%% Airfoil Numbers
% Values estimated from "Theory of Wing Sections" and "Fluid Dynamic Lift"

A0 = 0;             % Alpha where C_l = 0; [degree]
CL1_max = 1.4;      % max C_L pre-stall; [dimentionless]
CL2_max = 1.2;      % max C_L post-stall; [dimentionless]
S1 = 0.6/6;         % slope of C_l vs alpha graph at alpha = A0; [1/degrees]
CD0 = 0.005;        % CD at A0; [dimentionless]
CD1_max = 0.017;    % CD at ACD1; [dimentionless]
CD2_max = 1.8;      % CD at ACD2 (90 degrees); [dimentionless]
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
Chord = 1;                        %turbine blade chord; [m]
Planform_A = Blade_Length*Chord;    %planfrom blade area; [m^2]
AR = Blade_Length/Chord;            %aspect ratio; [dimentionless]
e = 0.70;

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


%% --------------- SIMULATE: ------------- %%

run_time = '10';
Run_time = str2double(run_time);

%running simulation
model = 'Three_blade_H_VAWT_Model';
load_system(model);
set_param(model, 'StopTime', run_time);
sim(model)


%% ------------ POST-PROCESSING: ------------- %%

animation_window_size = 25;
datum_aoa = 45; %degrees
datum = 100000;

% extract data:
blade_pos1_time = logsout.getElement('r_1').Values.Time;

blade_pos1 = logsout.getElement('r_1').Values.Data;
blade_pos2 = logsout.getElement('r_2').Values.Data;
blade_pos3 = logsout.getElement('r_3').Values.Data;

aoa1 = logsout.getElement('aoa 1').Values.Data;
aoa2 = logsout.getElement('aoa 2').Values.Data;
aoa3 = logsout.getElement('aoa 3').Values.Data;

L_I1 = logsout.getElement('L_I 1').Values.Data;
L_I2 = logsout.getElement('L_I 2').Values.Data;
L_I3 = logsout.getElement('L_I 3').Values.Data;
D_I1 = logsout.getElement('D_I 1').Values.Data;
D_I2 = logsout.getElement('D_I 2').Values.Data;
D_I3 = logsout.getElement('D_I 3').Values.Data;

blade1_torque = logsout.getElement('Torque1').Values.Data;
blade2_torque = logsout.getElement('Torque2').Values.Data;
blade3_torque = logsout.getElement('Torque3').Values.Data;
aero_torque_net = logsout.getElement('Aero_Torque_Net').Values.Data;
aero_torque_overlay = logsout.getElement('Torque_Overlay').Values.Data;

angular_velocity = logsout.getElement('Phi_dot').Values.Data;
angular_position = logsout.getElement('Phi').Values.Data;

blade1_force = logsout.getElement('F_net1').Values.Data;
blade2_force = logsout.getElement('F_net2').Values.Data;
blade3_force = logsout.getElement('F_net3').Values.Data;

% ------ plot ------ %
figure
subplot(2,1,1)
hold on
plot([0,Run_time],[0,0],'k-.')
plot(blade_pos1_time, aero_torque_net)
xlabel('Time [s]')
ylabel('Net Aerodynamic Torque [Nm]')
title('Net Aerodynamic Torque vs Time')

subplot(2,1,2)
hold on
plot([0,Run_time],[0,0],'k-.')
plot(blade_pos1_time, aero_torque_overlay)
xlabel('Time [s]')
ylabel('Aerodynamic Torque [Nm]')
title('Aerodynamic Torque Overlay of Three Blades')

figure
hold on
plot(blade_pos1_time, angular_velocity)
xlabel('Time [s]')
ylabel('Angular Velocity [rad/s]')
title('Angular Velocity of Turbine vs Time')

figure
subplot(2,1,1)
polarplot(angular_position, aoa1+datum_aoa,'b',angular_position, datum_aoa.*ones(length(angular_position)),'k-.',angular_position, (datum_aoa+ACD1).*ones(length(angular_position)),'r-.',angular_position, (datum_aoa-ACD1).*ones(length(angular_position)),'r-.')
title('Angle of Attack vs Phi')

subplot(2,1,2)
hold on
plot(blade_pos1_time, aoa1)
plot([0,Run_time],[0,0],'k-.')
plot([0,Run_time],[ACD1,ACD1],'r-.')
plot([0,Run_time],[-ACD1,-ACD1],'r-.')
xlabel('Time [s]')
ylabel('Angle of Attack [degrees]')
title('Single Blade Angle of Attack vs Time')

figure
polarplot(angular_position, aero_torque_net+datum,'b',angular_position, datum.*ones(length(angular_position)),'k-.')
title('Net Aero Torque vs Phi')

figure
subplot(1,2,1)
polarplot(angular_position, blade1_torque+datum/3,angular_position, datum/3.*ones(length(angular_position)),'k-.')
title('Single Blade Torque vs Phi')
subplot(1,2,2)
polarplot(angular_position, aoa1+datum_aoa,'b',angular_position, datum_aoa.*ones(length(angular_position)),'k-.',angular_position, (datum_aoa+ACD1).*ones(length(angular_position)),'r-.',angular_position, (datum_aoa-ACD1).*ones(length(angular_position)),'r-.')
title('Angle of Attack vs Phi')

figure
hold on
plot(blade_pos1_time, blade1_force + blade2_force + blade3_force)
xlabel('Time [s]')
ylabel('Force [N]')
title('Aero Force on Turbine vs Time')

%% Animation

%NACA 4-Digit Airfoil
%Airfoil
Airfoil = '0012';
chordlength = 1;
gridpoints = 20;

m  = str2double(Airfoil(1));
p  = str2double(Airfoil(2));
xx = str2double(Airfoil(3:4));

%define parameters
M = m/100;
P = p/10;
T = xx/100;
c = chordlength;
step = c/gridpoints;

%Mean Camber Line:
x = [0:step:c].';
y_c = [0];
dy_c = [0];
for j = 1:gridpoints + 1
   if x(j) >=0 && x(j) < P.*c;
        y_c(j,1)  = M/(P.^2).*x(j).*(2.*P - x(j)/c);
        dy_c(j,1) = 2.*M/(P.^2).*(P - x(j)/c);  
   else x(j) >= P.*c && x(j) <= c;
        y_c(j,1)  = M.*((c-x(j))/((1-P).^2)).*( 1+x(j)/c - 2.*P );
        dy_c(j,1) = 2.*M/((1 - P).^2).*(x(j)/c - P); 
   end
end

% figure
% hold on
% axis equal
% grid on
% plot(x,y_c,'b-');

%Thinkness Distribution
a_0 = 0.2969;
a_1 = -0.126;
a_2 = -0.3516;
a_3 = 0.2843;
a_4 = -0.1015;   %trailing edge pointy
%a_4 = -0.1036;   %trailing edge closed

y_t = T/0.2.*(  a_0.*x.^(0.5) + a_1.*x + a_2.*x.^2 + a_3.*x.^3 + a_4.*x.^4);

theta = atan(dy_c);
x_u = x - y_t.*sin(theta);
y_u = y_c + y_t.*cos(theta);
%plot(x_u,y_u,'k-')

x_L = x + y_t.*sin(theta);
y_L = y_c - y_t.*cos(theta);
%plot(x_L,y_L,'r-')
%title(Airfoil)

x_u = Chord.*x_u;
y_u = Chord.*y_u;
x_L = Chord.*x_L;
y_L = Chord.*y_L;

x_u = x_u - Chord/4;
x_L = x_L - Chord/4;

% foil_coordinates = [cat(1,x_u,flipud(x_L(2:end) ) ),(cat(1,y_u,flipud(y_L(2:end) ) ))];
% figure
% hold on
% plot(x_u,y_u,'b-');plot(x_L,y_L,'b-')
% axis equal

figure
for i = 1:length(blade_pos1_time)
   
    %set up window
    hold off
    plot(0,0,'kx');
    hold on; plot(animation_window_size,animation_window_size);
    plot(-animation_window_size,-animation_window_size);
    axis([-animation_window_size,animation_window_size,-animation_window_size,animation_window_size])
    grid on
    
    %---------- blade one ----------- %
    if blade1_torque(i) > 0
        blade1_draw_color = 'g';
    else
        blade1_draw_color = 'r';
    end
    %draw foil
    aoa_graphic_blade1 = -atan2(blade_pos1(i,2),blade_pos1(i,1)) + pi/2 + blade1_pitch.*pi/180;
    Upperfoil1 = [x_u,y_u]*Rotation_Matrix(aoa_graphic_blade1);
    x_u1 = Upperfoil1(:,1); y_u1 = Upperfoil1(:,2);
    Lowerfoil1 = [x_L,y_L]*Rotation_Matrix(aoa_graphic_blade1);
    x_L1 = Lowerfoil1(:,1); y_L1 = Lowerfoil1(:,2);
    
    plot(x_u1 + blade_pos1(i,1),y_u1 + blade_pos1(i,2),blade1_draw_color)
    plot(x_L1 + blade_pos1(i,1),y_L1 + blade_pos1(i,2),blade1_draw_color)
    
    %draw lift and drag
    plot([blade_pos1(i,1),blade_pos1(i,1)+L_I1(i,1)],[blade_pos1(i,2),blade_pos1(i,2)+L_I1(i,2)],'b-')
    plot([blade_pos1(i,1),blade_pos1(i,1)+D_I1(i,1)],[blade_pos1(i,2),blade_pos1(i,2)+D_I1(i,2)],'k-')
    
    
    %---------- blade two ----------- %
    if blade2_torque(i) > 0
        blade2_draw_color = 'g';
    else
        blade2_draw_color = 'r';
    end
    %draw foil
    aoa_graphic_blade2 = -atan2(blade_pos2(i,2),blade_pos2(i,1)) + pi/2 + blade2_pitch.*pi/180;
    Upperfoil2 = [x_u,y_u]*Rotation_Matrix(aoa_graphic_blade2);
    x_u2 = Upperfoil2(:,1); y_u2 = Upperfoil2(:,2);
    Lowerfoil2 = [x_L,y_L]*Rotation_Matrix(aoa_graphic_blade2);
    x_L2 = Lowerfoil2(:,1); y_L2 = Lowerfoil2(:,2);
    
    plot(x_u2 + blade_pos2(i,1),y_u2 + blade_pos2(i,2),blade2_draw_color)
    plot(x_L2 + blade_pos2(i,1),y_L2 + blade_pos2(i,2),blade2_draw_color)
    
    %draw lift and drag
    plot([blade_pos2(i,1),blade_pos2(i,1)+L_I2(i,1)],[blade_pos2(i,2),blade_pos2(i,2)+L_I2(i,2)],'b-')
    plot([blade_pos2(i,1),blade_pos2(i,1)+D_I2(i,1)],[blade_pos2(i,2),blade_pos2(i,2)+D_I2(i,2)],'k-')
    
    
    %---------- blade three ----------- %
    if blade3_torque(i) > 0
        blade3_draw_color = 'g';
    else
        blade3_draw_color = 'r';
    end
    %draw foil
    aoa_graphic_blade3 = -atan2(blade_pos3(i,2),blade_pos3(i,1)) + pi/2 + blade3_pitch.*pi/180;
    Upperfoil3 = [x_u,y_u]*Rotation_Matrix(aoa_graphic_blade3);
    x_u3 = Upperfoil3(:,1); y_u3 = Upperfoil3(:,2);
    Lowerfoil3 = [x_L,y_L]*Rotation_Matrix(aoa_graphic_blade3);
    x_L3 = Lowerfoil3(:,1); y_L3 = Lowerfoil3(:,2);
    
    plot(x_u3 + blade_pos3(i,1),y_u3 + blade_pos3(i,2),blade3_draw_color)
    plot(x_L3 + blade_pos3(i,1),y_L3 + blade_pos3(i,2),blade3_draw_color)
    
    %draw lift and drag
    plot([blade_pos3(i,1),blade_pos3(i,1)+L_I3(i,1)],[blade_pos3(i,2),blade_pos3(i,2)+L_I3(i,2)],'b-')
    plot([blade_pos3(i,1),blade_pos3(i,1)+D_I3(i,1)],[blade_pos3(i,2),blade_pos3(i,2)+D_I3(i,2)],'k-')
     
    drawnow
end