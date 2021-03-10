% Victoria (Vicky) Lopez 
% Function for ODE45 parameters

% Inputs: t->the time span of the function, initial_data-> xi, the
% initial conditions of the variables stated below.
% Outputs: the differential of the variables stated below throughout that
% time span.

function [dx] = ODE45_br(~,initial_data) 

   % Global variable definitions 
   global V_bottle  A_throat A_bottle V_air1 rho_air P_amb gamma R ...
       rho_water g CD Cd y0 l_s P_total M_air1 theta P_end
      
   % --------------------------------------------------------------- %
   
   %  xi(1) => initial_data(1) = x0, initial horizontal distance
   %  xi(2) => initial_data(2) = y0, initial vertical distance
   %  xi(3) => initial_data(3) = vx0, initial horizontal velocity
   %  xi(4) => initial_data(4) = vy0, initial vertical velocity
   %  xi(5) => initial_data(5) = V_air1, initial volume of air inside bottle
   %  xi(6) => initial_data(6) = M_air1, initial mass of the air inside the bottle
   %  xi(7) => initial_data(7) = M_rocket, mass of the bottle rocket
   
   % --------------------------- Phase 1 --------------------------- %
   
   % If statements to check if we are still on the stand
   heading = sqrt(initial_data(1)^2+(initial_data(2)-y0)^2);
   if heading < l_s
       hx = cos(theta);
       hz = sin(theta);
   else
       hx = initial_data(3)/sqrt(initial_data(3)^2+initial_data(4)^2);
       hz = initial_data(4)/sqrt(initial_data(3)^2+initial_data(4)^2);
   end
   
   % Compute pressure
     P = ((initial_data(6)/M_air1) ^ gamma) * P_end;

   if initial_data(5) < V_bottle
       
       % dvol_dt -> change in volume with time
       dx(5) = Cd * A_throat * sqrt((2 / rho_water) * (P_total * ...
           ((V_air1 / initial_data(5))^gamma ) - P_amb));
       
       % there is no change in the mass of the air
       dx(6) = 0;
       
       P = ((V_air1 / initial_data(5)) ^ gamma) * P_total;
       % dmassr_dt -> change in the mass of the rocket with time
       dx(7) = -Cd * A_throat * sqrt(2 * rho_water * (P - P_amb));
       
       % Compute thrust in phase 1
       T = 2 * Cd * A_throat * P_total;
       
    % --------------------------- Phase 2 --------------------------- %
    % initial_data(5) >= V_bottle && P > P_amb
    elseif  P > P_amb
        
       % there is no change in the volume 
       dx(5) = 0;
        
       rho = initial_data(6) / V_bottle;
       % Using the ideal gas law
       Temp = P/(rho*R);
       Pcritical = P * (2 / (gamma + 1))^(gamma / (gamma - 1));

       if Pcritical > P_amb

           T_e = (2 / (gamma + 1)) * Temp;
           P_e = Pcritical;
           V_e = sqrt(gamma * T_e * R);
           % Using the ideal gas law
           rho_e = Pcritical / (R * T_e);

       else
           
           Mach = sqrt(((P/P_amb)^(((gamma-1)/gamma))-1)*(2/(gamma-1)));
           %Mach = sqrt((nthroot(P/P_amb,gamma/(gamma-1))-1)/((gamma-1)/2));
           T_e = Temp * (1 + ((gamma - 1) / 2) * Mach ^2); 
           P_e = P_amb;
           V_e = Mach * sqrt(gamma * T_e * R);
           rho_e = P_amb / (R * T_e);


        end 
       
       flowrate = Cd * A_throat * V_e * rho_e ;
       
       % Compute thrust in phase 2
       T = flowrate * V_e + (P_e - P_amb) * A_throat;
       % dmassa_dt -> change in the mass of air with time
       dx(6) = -flowrate;
       % dmassr_dt -> change in the mass of the rocket with time
       dx(7) = flowrate;

   % --------------------------- Phase 3 --------------------------- %
   
       else
        
       % there is no change in the volume 
       dx(5) = 0;
       % there is no change in the mass of air with time
       dx(6) = 0;
       % there is no change in the mass of the rocket with time
       dx(7) = 0;
       % There is no thrust on phase 3
       T = 0;

   end
   
   V_total = sqrt(initial_data(3)^2+initial_data(4)^2);
   D = (rho_air / 2) * V_total^2 * CD * A_bottle;
   % Velocity 
   % dx_dt -> change in the horizontal distance
   dx(1) = initial_data(3);
   % dz_dt -> change in the vertical distance
   dx(2) = initial_data(4);

   % Acceleration in x and z
   % dvx_dt -> change in the horizontal velocity with time
   dx(3) = ((T - D) * hx) / initial_data(7);
   % dvz_dt -> change in the vertical velocity with time
   dx(4) = (((T - D) * hz)/initial_data(7)) - g;


   dx = dx';
   
end 



