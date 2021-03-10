clear all;
close all;
clc;

%% Define global variables
global g % acceleration due to gravity
global gamma % ratio of specific heats
global R % gas constant of air
global m_b % mass of empty 2-liter bottle with cones and fins
global v_b % volume of empty bottle
global c_d % discharge coefficient
global A_t % cross sectional area of throat
global p_w % density of water 
global p_a % ambient air density
global P_a % ambient/atmospheric pressure of air
global Pi_g % initial gage pressure of air in bottle
global Pi_0 % intial total pressure of air in rocket % One of 4 dependence parameters
global Ti_a % initial temperature of the air
global vi_w % initial volume of water inside bottle % One of the 4 dependence parameters
global vi_a % initial volume of air inside bottle
global P_end % the pressure of air in the bottle after the water has been expelled
global T_end % the temperature of air in the bottle after the water has been expelled
global mi_a % intial mass of air in the bottle
global C_D % coefficient of drag % One of 4 dependence parameters.
global A_B % cross section area of bottle
global x0 % initial x position of the rocket
global z0 % initial height of rocket, ie height of test stand
global l_rod % length of rod on test stand
global theta % initial angle of rocket on test stand
global Vx_0 % intial horizontal velocity of the rocket
global Vy_0 % initial vertical velocity of the rocket
global mi_R % initial mass of rocket

% Set g
g = 9.81; % m/s^2
% Set gamma
gamma = 1.4; % no units
% Set R
R = 287; % J/kgK
% Set m_b
m_b = 0.15; % kg
% Set v_b
v_b = 0.002; % m^3
% Set c_d
c_d = 0.8; % no units
% Set A_t
A_t = pi/4*(2.1/100)^2; % m^2
% Set p_w
p_w = 1000; % kg/m^3
% Set p_a
p_a = 0.961; % kg/m^3
% Set P_a
%P_a = psi2Pa(12.1);
P_a = 83426.56; % Pa
% Set P^i_g
%Pi_g = psi2Pa(50);
% Test case 1: Pi_g = 344738 + 10;
% Test case 2 : Pi_g = 344738 - 10;
% Test case 3: Pi_g = 344738 + 10000;
% Test case 4 : Pi_g = 344738 - 10000;
% Additional fine tuning was performed to get:
% Test case 5: Pi_g = 344738 + 21500;
Pi_g = 344738 + 21500; % Pa
% Set P^i_0
Pi_0 = P_a + Pi_g; % Pa
% Set T^i_a
Ti_a = 300; % K
% Set v^i_w
% Test case 1: vi_w = 0.0012;
% Test case 2: vi_w = 0.0015;
% Test case 3: vi_w = 0.0009;
% Test case 4: vi_w = 0.0006;
% Additional fine tuning was performed to get:
% Test case 5: vi_w = 0.00063; % Optimum
vi_w = 0.00063; % m^3
% Set v^i_a
vi_a = v_b-vi_w; % m^3
% Set P_end
P_end = Pi_0*((vi_a/v_b)^gamma); % Pa
% Set T_end
T_end = Ti_a*((vi_a/v_b)^(gamma-1)); % k
% Set mi_a
mi_a = (Pi_0*vi_a)/(R*Ti_a); % kg
% Test Case 1: C_D = 0.7;
% Test Case 2: C_D = 0.2;
% No additional fine tuning was conducted. However, it seems like it would
% be difficult to decrease the drag coefficient of the bottle rocket by
% more than 0.05. Thus I have set my drag coefficient to an arbitrary lower
% value of 0.45.
% Set C_D
C_D = 0.45;
% Set A_B
A_B = pi/4*(10.5/100)^2; % m^2
% Set x0
x0 = 0; % m
% Set z0
z0 = 0.25; % m 
% Set l_rod
l_rod = 0.5; % m
% Set theta
% Test Case 1: theta = 60;
% Test Case 2: theta = 90;
% Test Case 3: theta = 30;
% Test Case 4: theta = 20;
% No additional fine tuning required. 45 degrees is the already the optimum
% launch angle.
theta = 45; % degrees
theta = theta * (pi/180); % radians
% Set Vx_0
Vx_0 = 0; % m/s
% Set Vy_0
Vy_0 = 0; % m/s
% Set mi_R
mi_R = m_b+p_w*(v_b-vi_a)+(Pi_0/(R*Ti_a))*vi_a;

% % Make the global variables to be able to call them from every function 
% global V_bottle D_throat D_bottle M_bottle P_gage A_throat A_bottle V_air1...
%     rho_air P_amb gamma R T_air1 rho_water V_water1 g CD Cd v0 x0 y0 l_s ...
%     tspan vx0 vy0 P_total M_air1 theta M_water1 M_rocket P_end T_end
% 
% %Define the constants: 
% 
% % Bottle constants 
% V_bottle = 0.002;   % volume of empty bottle [m^3]
% D_throat = 0.021;   % diameter of throat [m] -> 2.1 [cm]
% D_bottle = 0.105;   % diameter of bottle [m] -> 10.5 [cm]
% M_bottle = 0.15;    % mass of the bottle [kg] 
% P_gage = 344738 + 21500;    % initial gage pressure of air in bottle [Pa] -> [psi]
% A_throat = (pi/4)*(D_throat)^2; % volume of throat [m^2]
% A_bottle = (pi/4)*(D_bottle)^2; % volume of bottle [m^2]
% 
% % Water constants 
% rho_water = 1000;   % density of water [kg/m^3]
% % V_water1 = 0.001;   % initial volume of water inside bottle [m^3]
% V_water1 = 0.00063;
% 
% % Ambient/atmospheric and air constants 
% rho_air = 0.961;    % ambient air density [kg/m^3]
% P_amb = 83426.56;   % atm pressure [Pa] -> 12.1 [psi]
% gamma = 1.4;    % ratio of specific heats for air
% R = 287;    % specific gas constant [J/kgK]
% T_air1 = 300;   % initial temperature of air [K]
% V_air1 = V_bottle - V_water1; % initial volume of air inside bottle [m^3]
% 
% % Force contants 
% g = 9.81;   % acceleration due to gravity [m/s^2]
% CD = 0.45;   % drag coefficient
% Cd = 0.8;   % discharge coefficient
% P_total = P_amb + P_gage; % Total pressure [Pa]
% 
% P_end = P_total*((V_air1/V_bottle)^gamma);
% 
% T_end = T_air1*((V_air1/V_bottle)^(gamma-1));
% 
% % Computation constants and initial conditions 
% v0 = 0.0;    % initial velocity of rocket [m/s]
% % split V into vy0 and vx0
% vx0 = 0;
% vy0 = 0;
% x0 = 0.0;    % initial horizontal distance [m]
% y0 = 0.25;  % initial vertical height [m]
% l_s = 0.5;  % length of test stand [m]
% tspan = [0 5];  % integration time [s]
% 
% % max_h = 17.37;   % maximum height [m]
% % max_dist = 60.45;      % maximum distance [m]
% 
% % ------ Test cases for theta ------ %
% % Test 1: 30
% % Test 2: 60
% % Test 3: 90
% % Optimum 45º
% theta = 45;  % Initial angle of rocket 
% % change it to radians
% theta = theta*(pi/180);
% 
% % compute the initial mass of the air inside the bottle
% M_air1 = (P_total*V_air1)/(R*T_air1);
% % compute the initial mass of the water inside the bottle
% M_water1 = rho_water * V_water1;
% % compute the initial mass of the bottle rocket 
% M_rocket = M_bottle + rho_water * ((V_bottle - V_air1) + (P_total / R * T_air1)) * V_air1;
% 





% Build initial conditions vector
xi = [x0,z0,Vx_0,Vy_0,vi_a,mi_a,mi_R];

% Call the rocket function to analyze the trajectory of the bottle rocket
[t,x] = ode45('rocket',[0 5],xi);

% Parse output vector for ease of graphing
d = x(:,1); % distance (x position)
h = x(:,2); % height (z position)
Vx = x(:,3); % velocity (x component)
Vz = x(:,4); % velocity (z component)
v_a = x(:,5); % volume of air
m_a = x(:,6); % mass of air
m_R = x(:,7); % mass of rocket

% Call ThrustVec function to produce set of thrust values for plotting
% F = ThrustVec(v_a,m_a);
% 
% find the maximum distance
maxd = max(d(h>0)); % restricts maximum height and distance based on the criteria that the height cannot be < 0
% find the maximum height
maxh = max(h(h>0));

% Plot Height vs. Distance
figure(1)
plot(d(h>-1),h(h>-1),'LineWidth',2) % d(h>-1) was chosen to allow the trajectory to cross the x-axis. ODE45 produced a point at location x where y is slightly > 0 and one where y is slightly < 0, but no exact x value exists for y = 0 
title('Height vs Distance');
ylabel('Height (m)');
xlabel('Distance (m)');
