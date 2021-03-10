% Andrew Logue
% 12/6/2020
% Project 2 - ASEN 2012

clc
clear
close all

% initialize global constants - for use in external functions
global g Cd airDensity emptyBottleVolume atmosphericPressure... 
    specificHeatRatio waterDensity throatDiameter bottleDiameter R ...
    bottleMass CD gagePressure waterVolume T_0 vx_0 vz_0 theta x_0 z_0...
    throatArea bottleArea totalPressure totalRocketMass_0 tspan...
    waterMass airMass_0 airVolume_0 standL T_f P_f
    

% define global constants
g = 9.81; % m/s2 acceleration due to gravity
Cd = 0.8; % discharge coefficient
airDensity = 0.961; % kg/m^3 ambient air density
waterDensity = 1000; % kg/m^3  density of water
atmosphericPressure = 12.1 * 6894.76; % Pa ambient atmospheric pressure
gagePressure = 377000; % Pa  initial gage pressure of air in bottle
totalPressure = atmosphericPressure + gagePressure; % Pa total pressure
R = 287; % J/kgK  gas constant of air
CD = 0.3; %  drag coefficient
T_0 = 300; % K  initial temperature of air
tspan = [0:0.001:5]; % integration time = 0sec to 5sec

% bottle properties
emptyBottleVolume = 0.002; % m^3  volume of empty bottle
specificHeatRatio = 1.4; % ratio of specific heats for air
throatDiameter = .021; % m  diameter of throat
throatArea = (pi/4) * throatDiameter^2; % m^2 area of throat
bottleDiameter = .105; % m  diameter of bottle
bottleArea = (pi/4) * bottleDiameter^2; % m^2 area of bottle
bottleMass = 0.15; % kg  mass of empty 2-liter bottle with cone and fins
waterVolume = 0.001; % m^3  initial volume of water inside bottle
airVolume_0 = emptyBottleVolume -...
    waterVolume; % m^3  initial volume of air inside bottle
waterMass = waterDensity*waterVolume; % kg mass of water inside the bottle
airMass_0 = (totalPressure * airVolume_0) /...
    (R * T_0); % kg mass of air inside the bottle
totalRocketMass_0 = bottleMass +...
    waterMass + airMass_0; % kg total mass of the bottle rocket

% initial launch conditions
vx_0 = 0.0; % m/s  initial velocity of rocket (vertical)
vz_0 = 0.0; % m/s  initial velocity of rocket (horizontal)
theta = 50 * pi / 180; % radians  initial angle of rocket
x_0 = 0.0; % m  initial horizontal distance
z_0 = 0.25; % m  initial vertical height
standL = 0.5; % m length of test stand

% final launch conditions
% K final Temperature
T_f = T_0 * ((airVolume_0/emptyBottleVolume)^(specificHeatRatio-1));
% Pa final pressure
P_f = totalPressure * ((airVolume_0/emptyBottleVolume)^specificHeatRatio);

% initial conditions vector
x_i = [x_0, z_0, vx_0, vz_0, airVolume_0, airMass_0, totalRocketMass_0];
% ode45 function
[t, x] = ode45('fun', tspan, x_i);

% parse data for graphing 
x_f = x(:,1); % distance (x position)
z_f = x(:,2); % height (z position)
vx_f = x(:,3); % velocity (x component)
vz_f = x(:,4); % velocity (z component)
airVolume_f = x(:,5); % volume of air
airMass_f = x(:,6); % mass of air
totalRocketMass_f = x(:,7); % mass of rocket

% find thrust
[thrust_1 ,thrust_2, thrust_3] = thrust(x);

% find the maximum height (where height must be greater than 0)
max_z_f = max(z_f(z_f > 0));
% find the maximum distance (where height must be greater than 0)
max_x_f = max(x_f(z_f > 0));

figure(1)
plot(x_f(z_f > -1), z_f(z_f > -1), 'LineWidth', 1.5)
xlim([0 (max_x_f + 5)])
ylim ([0 (max_z_f + 5)])
grid on
title('Bottle Rocket Travel: Height vs Distance');
ylabel('Height (meters)');
xlabel('Distance (meters)');

% % second plot, the thrust of the botttle rocket vs. time 
figure(2);
plot(t,[thrust_1 thrust_2 thrust_3],'LineWidth',2) 
hold on
grid on
xlim([0 1/2])
title('Bottle Rocket Thrust vs Time')
xlabel('Time (Seconds)')
ylabel('Thrust (Newtons)')

% ADDITIONAL VERIFICATION PLOTS
figure(3)
plot(t, vx_f, 'LineWidth', 1.5)
title('velocity (x component) vs Time')
xlabel('Time (seconds)')
ylabel('Velocity in xndirection (meters per second)')

figure(4)
plot(t, vz_f, 'LineWidth', 1.5)
title('velocity (z component) vs Time')
xlabel('Time (seconds)')
ylabel('Velocity in z direction (meters per second)')

figure(5)
plot(t, airVolume_f, 'LineWidth', 1.5)
title('Change in Bottle Rocket Volume of Air')
xlim([0 0.25])
ylim([0.001 0.002])
xlabel('Time (Seconds)')
ylabel('Volume of Air (meters^3)')
