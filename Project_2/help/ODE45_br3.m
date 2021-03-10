% Victoria (Vicky) Lopez 
% Bottle rocket function used in ode45 for Project 2
% Date created: Nov 12th, 2020
% Last modified on: Dec 1st, 2020

% Inputs: - x ==> initial conditions vectors (xi vector from main function)
%           later explained in detail under 'Initial conditions'.
%         - t ==> variable for time (it is not used in the function, could
%         be replaced with '~' but wasn't for consistency)

% ------------------------ Initial conditions ------------------------ %
   
%  xi(1) => initial_data(1) = x0, initial horizontal distance
%  xi(2) => initial_data(2) = y0, initial vertical distance
%  xi(3) => initial_data(3) = vx0, initial horizontal velocity
%  xi(4) => initial_data(4) = vy0, initial vertical velocity
%  xi(5) => initial_data(5) = V_air1, initial volume of air inside bottle
%  xi(6) => initial_data(6) = M_air1, initial mass of the air inside the bottle
%  xi(7) => initial_data(7) = M_rocket, mass of the bottle rocket

% -------------------------------------------------------------------- % 

% Outputs: - dx ==> derivative vector - later explained in detail under 
%            'Derivatives assignment'.
            
% --------------------- Derivatives assignment ----------------------- %
   
%  dx(1) => change in the horizontal distance (horizontal velocity)
%  dx(2) => change in the vertical distance (vertical velocity)
%  dx(3) => change in the horizontal velocity (horizontal acceleration)
%  dx(4) => change in the vertical velocity (vertical acceleration)
%  dx(5) => change in the volume of air inside bottle
%  dx(6) => change in the mass of the air inside the bottle
%  dx(7) => change in the mass of the bottle rocket

% -------------------------------------------------------------------- % 

% Function explication: this functin will compute the change in the
% conditions stated before, meaning it will take the derivative of each
% equation. 


function [dx] = ODE45_br3(t,initial_data) 

   % Global variable definitions 
   global V_bottle  A_throat A_bottle V_air1 rho_air P_amb gamma R...
       rho_water g CD Cd y0 l_s M_air1 theta P_gage P_end
   
   % ----------------------------- Phase 1 ------------------------------ %
   
   % Compute the total pressure (gage pressure plus the ambient pressure)
   % since it's a value that will be used frequently throughout the function
   P_total = P_amb + P_gage;
   
   % Check if the rocket started the first phase (the initial volume of air
   % inside bottle is less than the volume of the bottle)
   if initial_data(5) < V_bottle
       
        % Check if we still on the test stand (l_s):
        
        h = sqrt((initial_data(1)^2)+(initial_data(2)-y0)^2);

        if h <= l_s
            % using Pythagorean Theorem  
            v_total = sqrt((initial_data(4).^2) + (initial_data(3).^2));
            hx = cos(theta);
            hz = sin(theta);
        else
            v_total = sqrt((initial_data(4).^2) + (initial_data(3).^2));
            hx = initial_data(3)/v_total;
            hz = initial_data(4)/v_total;
        end
        
        
        % Get the 3 forces acting on the rocket in phase 1: 
        % DRAG 
        D = (rho_air / 2) .* (v_total).^2 * CD * A_bottle;
        % PRESSURE
        P = (( V_air1 ./ initial_data(5)) .^ gamma) .* P_total;
        % THRUST
        Th = 2.* Cd .* A_throat .* (P - P_amb);
        
        
   % ------------------------- derivatives P1 --------------------------- %
        
        % Velocity 
        % dx(1) -> change in the horizontal distance
        dx_dt = initial_data(3);
        % dx(2) -> change in the vertical distance
        dz_dt = initial_data(4);

        % Acceleration in x and z
        % dx(3) -> change in the horizontal velocity with time
        dvx_dt = ((Th - D) * hx) ./ initial_data(7);
        % dx(4) -> change in the vertical velocity with time
        dvz_dt = (((Th - D) * hz) - initial_data(7) * g)./ initial_data(7);

        % dx(5) -> change in volume with time
        dvol_dt = Cd * A_throat * sqrt((2 / rho_water) * (P_total * ...
           ((V_air1 / initial_data(5))^gamma ) - P_amb));

        % dx(6) -> there is no change in the mass of the air
        dmassa_dt = 0;

        % dx(7) -> change in the mass of the rocket with time
        dmassr_dt = -Cd .* A_throat .* sqrt(2 .* rho_water .* (P - P_amb));
            
        % Make the output a column vector 
        dx = [dx_dt;dz_dt;dvx_dt;dvz_dt;dvol_dt;dmassa_dt;dmassr_dt];

   % ----------------------------- Phase 2 ------------------------------ %
   
   % Check if the rocket started the second phase (the volume of air
   % inside bottle is equal to the volume of the bottle) which means the
   % bottle ran out water.
   
   elseif initial_data(5) >= V_bottle
        
       % Definitely the rocket isn't on the test stand, therefore:
        v_total = sqrt((initial_data(4).^2) + (initial_data(3).^2));
        hx = initial_data(3)/v_total;
        hz = initial_data(4)/v_total;
        
        P_state = P_end * (initial_data(6)/M_air1)^(gamma);

        % Check condition that separates the 2nd phase from the 3rd
        % (preassure) 
         if P_state > P_amb

            P_cr = (P_state) * (2./(gamma+1)).^(gamma/(gamma-1));
            rho = initial_data(6) / V_bottle;
            T = P_state/(rho*R);
            
  % Compute the exit values for each condition using the equations from the equation sheet
            if P_cr <= P_amb
                
               Mach = sqrt(((P_state/P_amb)^(((gamma-1)/gamma))-1)*(2/(gamma-1)));
               T_e = T/(1+((gamma-1)/2)*Mach^2);
               P_e = P_amb;
               rho_e = P_amb/(R*T_e) ;
               V_e = Mach * sqrt(gamma*T_e*R);
               
            elseif P_cr > P_amb
                
                P_e = P_cr;
                T_e = (2/(gamma+1))*T;
                V_e = sqrt(gamma*T_e*R);
                rho_e = P_cr/(R*T_e);

             end

        % Compute the air's mass flow rate
        flowrate = Cd * rho_e * A_throat * V_e;
  
        % Get the 3 forces acting on the rocket in phase 2:
        % THRUST
        Th = flowrate *V_e + (P_e - P_amb) * A_throat;
        % DRAG
        D = (rho_air / 2) .* (v_total).^2 * CD * A_bottle; 

   % -------------------------- Derivatives P2 -------------------------- %

       % Velocity 
       % dx(1) -> change in the horizontal distance
       dx_dt = initial_data(3);
       % dx(2) -> change in the vertical distance
       dz_dt = initial_data(4);

       % Acceleration in x and z
       % dx(3) -> change in the horizontal velocity with time
       dvx_dt = ((Th - D) * hx) ./ initial_data(7);
       % dx(4) -> change in the vertical velocity with time
       dvz_dt = (((Th - D) * hz) - initial_data(7)*g) ./ initial_data(7);
       
       % dx(5) -> there is no change in the volume 
       dvol_dt = 0;
       % dx(6) -> change in the mass of air with time
       dmassa_dt = -flowrate;
       % dx(7) -> change in the mass of the rocket with time
       dmassr_dt = -flowrate;
        
       % Return the derivatives as a column vector for the ode call
       dx = [dx_dt;dz_dt;dvx_dt;dvz_dt;dvol_dt;dmassa_dt;dmassr_dt];
         
  % ----------------------------- Phase 3 ------------------------------ %
        
  % if the pressure inside the rocket is greater than the ambient
  % pressure the the rocket entered ballistic phase. 
  
        else

            % Definitely the rocket isn't on the test stand, therefore: 
            v_total = sqrt( (initial_data(4).^2) + (initial_data(3).^2) );
            hx = initial_data(3)/v_total;
            hz = initial_data(4)/v_total;
            
            % Get the forces acting on the rocket in phase 3:
            % THRUST
            Th = 0 ; % There is no thurst in the ballistic phase 
            % DRAG
            D = ( rho_air / 2) .* (v_total).^2 * CD*A_bottle; 

 % --------------------------- Derivatives P3 --------------------------- %
 
       % Velocity 
       %dx(1) -> change in the horizontal distance
       dx_dt = initial_data(3);
       % dx(2) -> change in the vertical distance
       dz_dt = initial_data(4);

       % Acceleration in x and z
       % dx(3) -> change in the horizontal velocity with time
       dvx_dt = ((Th - D) * hx) ./ initial_data(7);
       % dx(4) -> change in the vertical velocity with time
       dvz_dt = (((Th - D) * hz) - initial_data(7)*g) ./ initial_data(7);
       
       % dx(5) -> there is no change in the volume with time 
       dvol_dt = 0;
       % dx(6) -> there is no change change in the mass of air with time
       dmassa_dt = 0;
       % dx(7) -> there is no change change in the mass of the rocket with time
       dmassr_dt = 0;
        
       % Return the derivatives as a column vector for the ode call
       dx = [dx_dt;dz_dt;dvx_dt;dvz_dt;dvol_dt;dmassa_dt;dmassr_dt];

        end   
   end
   
end 