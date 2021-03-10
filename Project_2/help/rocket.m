function [dx] = rocket(~,x)
% x(1) = x, position in x direction
% x(2) = z, position in z direction
% x(3) = vx, velocity in x direction
% x(4) = vz, velocity in z direction
% x(5) = vair, volume of air in the rocket
% x(6) = mair, mass of air in the rocket
% x(7) = mR, mass of rocket

% dx(1) = Vx, change in position in x direction with time
% dx(2) = Vz, change in position in z direction with time
% dx(3) = ax, acceleration in x direction 
% dx(4) = az, acceleration in z direction
% dx(5) = dvairdt, change in volume of air in the rocket with time
% dx(6) = dmairdt, change in mass of air in the rocket with time
% dx(7) = dmRdt, change in mass of rocket with time

% Import all necessary global variables
global g % acceleration due to gravity
global gamma % ratio of specific heats
global R % gas constant of air
global v_b % volume of empty bottle
global c_d % discharge coefficient
global A_t % cross sectional area of throat
global p_w % density of water
global p_a % ambient air density
global P_a % ambient/atmospheric pressure of air
global Pi_0 % intial total pressure of air in rocket
global vi_a % initial volume of air inside bottle
global P_end % the pressure of air in the bottle after the water has been expelled
global mi_a % intial mass of air in the bottle
global C_D % coefficient of drag
global A_B % cross sectional area of bottle
global z0 % initial height of rocket, ie height of test stand
global l_rod % length of rod on test stand
global theta % initial angle of rocket on test stand

% Determine the heading of the rocket
if sqrt(x(1)^2+(x(2)-z0)^2) < l_rod
    % On test stand
    hx = cos(theta); % because this is the initial angle of the rocket
    hz = sin(theta);
else
    % Not on test stand
    hx = x(3)/sqrt(x(3)^2+x(4)^2);
    hz = x(4)/sqrt(x(3)^2+x(4)^2);
end

% Calculate pressure of air inside rocket assuming thrust stage 2 or
% ballistic phase. This pressure will be compared to P_a to determine
% whether in thrust stage 2 or ballistic phase.
P = ((x(6)./mi_a)^gamma)*P_end;

% Determine if in thrust phase 1
if x(5) < v_b
    % Set dvairdt using equation 10
    dx(5) = c_d*A_t*sqrt((2/p_w)*(Pi_0*((vi_a/x(5))^gamma)-P_a));
    % Set dmairdt
    dx(6) = 0; % Mass of air in bottle doesn't change in phase 1
    % Set dmRdt using equation 3 and 11
    P = ((vi_a/x(5))^gamma)*Pi_0;
    dx(7) = -c_d*A_t*sqrt(2*p_w*(P-P_a));
    % Find the thrust produced by the rocket in phase 1
    F = 2*c_d*(P-P_a)*A_t;
% Determine if in thrust phase 2
elseif P > P_a
    % Set dvairdt
    dx(5) = 0; % Once all of the water has been expelled, the volume of air inside the bottle will not change.
    % Find density of air inside the bottle
    p = x(6)/v_b;
    % Find temperature of air inside the bottle
    T = P/(p*R);
    % Find critical pressure
    P_star = P*((2/(gamma+1))^(gamma/(gamma-1)));
    % Determine if flow is choked
    if P_star > P_a
        % Find the exit pressure
        P_e = P_star;
        % Find the exit temperature
        T_e = (2/(gamma+1))*T;
        % Find the exit density
        p_e = P_e/(R*T_e);
        % Find exit velocity
        V_e = sqrt(gamma*R*T_e);
    else
        % The flow is not choked
        P_e = P_a;
        % Find the exit Mach number
        M_e = sqrt((nthroot(P/P_a,gamma/(gamma-1))-1)/((gamma-1)/2));
        % Find the exit temperature
        T_e = T*(1+((gamma-1)/2)*M_e^2);
        % Find the exit density
        p_e = P_a/(R*T_e);
        % Find the exit velocity 
        V_e = M_e*sqrt(gamma*R*T_e);
    end
    % Find dmairdt, the rate of change of the mass of air in the rocket
    dx(6) = c_d*p_e*A_t*V_e;
    % Find the thrust of the bottle rocket
    F = dx(6)*V_e+(P_e-P_a)*A_t;
    % Change the sign of dmairdt so that it is negative. Note that dmairdt
    % has to be positive when used to calculate thrust in order to produce
    % a positive thrust (negative thrust is meaningless). On the other
    % hand, we known from intuition that dmairdt should actually be
    % NEGATIVE, because the mass of air in the bottle is decreasing as its
    % being jettisoned out the end as thrust. Hence the sign change. 
    dx(6) = -dx(6);
    % Find dmRdt, the rate of change of the mass of the rocket
    % By intuition, dmRdt should also be negative. This is because air is 
    % being jettisoned FROM the rocket, which is thus losing mass. If it is
    % losing mass, then dmRdt is negative.
    dx(7) = dx(6);
else
    % Rocket is in ballistic phase
    % Set dvairdt
    dx(5) = 0; % Once all of the water has been expelled, the volume of air inside the bottle will not change.
    dx(6) = 0; % In the ballistic phase, air is no longer being expelled from the rocket
    dx(7) = 0; % If neither air nor water is being expelled from the rocket, its mass will not change.
    F = 0; % In the ballistic phase, no thrust is being produced. Hence the name.
end

% Calculate the drag acting on the rocket
V_R = sqrt(x(3)^2+x(4)^2);
D = (p_a/2)*V_R^2*C_D*A_B;
% Define dxdt = vx
dx(1) = x(3);
% Define dzdt = vz
dx(2) = x(4);
% Define ax
dx(3) = ((F-D)*hx)/x(7);
% Define az
dx(4) = (((F-D)*hz)/x(7))-g;

% ode45 requires output of a column vector
dx = dx';
end