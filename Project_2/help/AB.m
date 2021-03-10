
clear;
clc;
close all;


%% inital conditions 

global TestStandLength Theta Pgage Pamb Cd ThroatArea CD BottleArea Rhoairamb...
    RhoWater Volbottle y0 VAirInit GammaGas g TAirInit MassAirInit R

g = 9.81; % m/s2, acceleration due to gravity,
Cd= 0.8; % discharge coefficient
Rhoairamb = 0.961; % kg/m^3 ambient air density
Volbottle= 0.002; % m^3 volume of empty bottle
Pamb= 12.1*6894.76; % converted to Pa, atmospheric pressure
GammaGas = 1.4; % ratio of specific heats for air
RhoWater = 1000; % kg/m^3, density of water
DThroat= 2.1; % cm, diameter of throat
DBottle= 10.5; % in cm, diameter of bottle
R = 287; %J/kgK, gas constant of air
MBottle= 0.15; % kg mass of empty 2-liter bottle with cone and fins
CD= 0.5; % drag coefficient
Pgage= 50*6894.76; % in pascal, the 6894.76 is to convert. initial gage pressure of air in bottleVolwater,
VWaterInit= 0.001; % m^3, initial volume of water inside bottle
TAirInit = 300; % K, initial temperature of
Airv0 = 0.0 ;% m/s, initial velocity of rocket
Theta= 45 ; % initial angle of rocket in degress
X0 = 0.0; % in meters, initial horizontal distance
y0 = 0.25; % in m, initial vertical height
TestStandLength= 0.5; % in m, length of test stand
VAirInit = Volbottle - VWaterInit ; %initial volume of Air.
ThroatArea = pi * ((DThroat*10^-2)/2)^2; %Area of throat
BottleArea  = pi * ((DBottle*10^-2)/2)^2; %Bottle Area
TotalMass0 = MBottle + (VWaterInit*RhoWater) + (((Pgage+Pamb)*VAirInit ) / (R*TAirInit)); % Total mass
MassAirInit = (((Pgage+Pamb)*VAirInit ) / (R*TAirInit));

global t2 t1 t3
%globals are where the time of each phase ends, so last elment in t1 is
%where phase 1 ends, so on.

%% assign golbals:

t1 = [0 0];
t2 = [0];
t3 = [0];
%% 

%initial conditions:

VelX0 = 0;
VelZ0 = 0;
Range0 = 0;
Height0 = y0;

tspan = [0 5];

% xi = [TotalMass0 MassAirInit VAirInit VelX0 VelZ0 Range0 y0 ]; 
xi = [ Range0 y0 VelX0 VelZ0 VAirInit  MassAirInit TotalMass0];
% Call ODE
[ Time,  Results ] = ode45(@(Time,States) ODE45_br3(Time,States),tspan,xi);

d = Results(:,1); % distance (x position)
h = Results(:,2); % height (z position)
Vx = Results(:,3); % velocity (x component)
Vz = Results(:,4); % velocity (z component)
v_a = Results(:,5); % volume of air
m_a = Results(:,6); % mass of air
m_R = Results(:,7); % mass of rocket

% Call ThrustVec function to produce set of thrust values for plotting
% F = ThrustVec(v_a,m_a);

% % % Only to check with the test case
% % find the maximum distance
% maxd = max(d(h>0)); % restricts maximum height and distance based on the criteria that the height cannot be < 0
% % find the maximum height
% maxh = max(h(h>0));

% Plot Height vs. Distance
figure(1)
plot(d(h>-1),h(h>-1),'LineWidth',2,'Color',[0.3 0.7 1]) % d(h>-1) was chosen to allow the trajectory to cross the x-axis. ODE45 produced a point at location x where y is slightly > 0 and one where y is slightly < 0, but no exact x value exists for y = 0 
xlim([0 65])
title('Height vs Distance');
ylabel('Height (m)');
xlabel('Distance (m)');



%% calculate thrust

%ODE Has really a hard way of calculating the time, even though matrices t1
%and t2 can tell us the end time of each phase there's no time that matches
%it exactly in the Time matrix, and extracting the force of thrust via
%globals really yields bad results, thus the forces has to be calculated
%manually.


%index namings: those will use to index the storing of thrust and time
store1 = 0;
store2 = 0;
store3 = 0;

[ r , c ] = size(Results);
for i=1:r
    
if Results(i,5) < Volbottle
    
Pressure1(i) = ( ( VAirInit ./ Results(i,5) ) .^ GammaGas ) .* (Pgage+Pamb) ; 
Thrust1(i) = 2.* Cd .* ThroatArea .* ( Pressure1(i) - Pamb) ;
TP1(i) = Time(i);

%phase 2
elseif Results(i,5)>= Volbottle
    
    %T and P of end states
Tend = TAirInit * (( VAirInit/Volbottle) ^ (GammaGas-1) );
Pend = (Pgage+Pamb) * (( VAirInit/Volbottle) ^ (GammaGas) );
    
PressureCond = Pend * (Results(i,6)/MassAirInit)^(GammaGas) ;

if PressureCond>Pamb
Density = Results(i,6) / Volbottle;
Temp = PressureCond/(Density*R);
CriticalP = (PressureCond) * (2./(GammaGas+1)).^(GammaGas/(GammaGas-1));
  
if CriticalP > Pamb
    
    Mach  = 1;
    Texit = (2/(GammaGas+1))*Temp ;
    Vexit = sqrt(GammaGas*Texit*R);
    Pexit = CriticalP;
    Densityexit = CriticalP/(R*Texit) ;
    
elseif CriticalP <= Pamb
    
   Mach = sqrt(( (PressureCond/Pamb)^( ( (GammaGas-1)/GammaGas)) - 1 ) * (2/(GammaGas-1)));
   Texit = Temp/(1+((GammaGas-1)/2)*Mach^2);
   Pexit = Pamb;
   Densityexit = Pamb/(R*Texit) ;
   Vexit = Mach * sqrt(GammaGas*Texit*R);
end

% how mass of air and rocket change with time

MassAirFlowRate = Cd*Densityexit*ThroatArea*Vexit;

store2 = store2 + 1;
Thrust2(store2) = MassAirFlowRate *Vexit + (Pexit-Pamb)*ThroatArea ;
TP2(store2) = Time(i);

else
store3 = store3 + 1;
Thrust3(store3) = 0 ;
TP3(store3) = Time(i);


%% Phase 3: 
end

% NO THRUST



end

end


%% plot thrust:

figure(2);
plot([TP1 TP2 TP3],[Thrust1 Thrust2 Thrust3],'Color',[0 0 0],'LineWidth',2)
hold on
plot(TP1(end),Thrust1(end),'*','Color',[1 0 0],'MarkerSize',10,'MarkerFaceColor',[1 0 0])
plot(TP2(end),Thrust2(end),'*','Color',[0 0 1],'MarkerSize',10,'MarkerFaceColor',[0 0 1])
xlim([0 1/2])
grid minor
title('Thrust VS Time')
xlabel('Time (Seconds)')
ylabel('Thrust (N)')
legend('Thrust','End of Water Phase','End of Air Phase','Location','NorthEast')

