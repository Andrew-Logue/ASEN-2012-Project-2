function [dx] = ODE45_br2(t,initial_data) 

   % Global variable definitions 
   global V_bottle  A_throat A_bottle V_air1 rho_air P_amb gamma R...
       rho_water g CD Cd y0 l_s M_air1 theta t1 t2 t3 P_gage 

   % --------------------------------------------------------------- %
   
   %  xi(1) => initial_data(1) = x0, initial horizontal distance
   %  xi(2) => initial_data(2) = y0, initial vertical distance
   %  xi(3) => initial_data(3) = vx0, initial horizontal velocity
   %  xi(4) => initial_data(4) = vy0, initial vertical velocity
   %  xi(5) => initial_data(5) = V_air1, initial volume of air inside bottle
   %  xi(6) => initial_data(6) = M_air1, initial mass of the air inside the bottle
   %  xi(7) => initial_data(7) = M_rocket, mass of the bottle rocket
   
   % --------------------------- Phase 1 --------------------------- %
   
   P_total = P_gage + P_amb;
   % Compute pressure
   if initial_data(5) < V_bottle
       
       % If statements to check if we are still on the stand
       heading = sqrt(initial_data(1)^2+(initial_data(2)-y0)^2);
       if heading <= l_s
           % Compute lengths with pythagorean theorem
           V_total = sqrt((initial_data(4).^2) + (initial_data(3).^2));
           hx = cos(theta);
           hz = sin(theta);
       else
           V_total = sqrt((initial_data(3).^2)+(initial_data(4).^2));
           hx = initial_data(3)/sqrt(initial_data(3)^2+initial_data(4)^2);
           hz = initial_data(4)/sqrt(initial_data(3)^2+initial_data(4)^2);
       end

       
       % Get the 3 forces acting on the rocket in phase 1:
       
       % Compute drag
       D = (rho_air / 2) .* (V_total).^2 * CD * A_bottle;
       
       % Compute pressure
       P = ((V_air1 ./ initial_data(5)) .^ gamma) .* P_total;
       
       % Compute thrust 
       Th = 2 .* Cd .* A_throat .* (P - P_amb);
        
       % Velocity 
       % dx_dt -> change in the horizontal distance
       dx_dt = initial_data(3);
       % dz_dt -> change in the vertical distance
       dz_dt = initial_data(4);
       
       % Acceleration in x and z
       % dvx_dt -> change in the horizontal velocity with time
       dvx_dt = ((Th - D) * hx) ./ initial_data(7);
       % dvz_dt -> change in the vertical velocity with time
       dvz_dt = (((Th - D) * hz) - initial_data(7) * g)./ initial_data(7);
       
       % dvol_dt -> change in volume with time
       dvol_dt = Cd * A_throat * sqrt((2 / rho_water) * (P_total * ...
           ((V_air1 / initial_data(5))^gamma ) - P_amb));
       
       % there is no change in the mass of the air
       dmassa_dt = 0;
       
       % dmassr_dt -> change in the mass of the rocket with time
       dmassr_dt = -Cd .* A_throat .* sqrt(2 .* rho_water .* (P - P_amb));
       
       dx = [dx_dt;dz_dt;dvx_dt;dvz_dt;dvol_dt;dmassa_dt;dmassr_dt];
       t1 = [t1 ; [t Th]];
       
    % --------------------------- Phase 2 --------------------------- %
    % initial_data(5) >= V_bottle && P > P_amb
    elseif  initial_data(5) >= V_bottle 
        
%         % final temperature
%        T_f = T_air1 * ((V_air1/V_bottle)^(gamma-1));
       
       % final pressure
       P_f = P_total * ((V_air1/V_bottle)^(gamma));
        
       V_total = sqrt((initial_data(3).^2)+(initial_data(4).^2));
       hx = initial_data(3)/V_total;
       hz = initial_data(4)/V_total;
       
       
       P_s = P_f * (initial_data(6)/M_air1)^gamma;
       
       if P_s > P_amb
           
           rho = initial_data(6) / V_bottle;
           % Using the ideal gas law
           Temp = P_s/(rho*R);
           Pcritical = P_s * (2 ./ (gamma + 1)).^(gamma / (gamma - 1));
           
           if Pcritical > P_amb

               T_e = (2 / (gamma + 1)) * Temp;
               P_e = Pcritical;
               V_e = sqrt(gamma * T_e * R);
               % Using the ideal gas law
               rho_e = Pcritical / (R * T_e);

           elseif Pcritical <= P_amb
               
               Mach = sqrt(((P_s/P_amb)^(((gamma-1)/gamma))-1)*(2/(gamma-1)));
               %Mach = sqrt((nthroot(P/P_amb,gamma/(gamma-1))-1)/((gamma-1)/2));
               T_e = Temp / (1 + ((gamma - 1) / 2) * Mach ^2); 
               P_e = P_amb;
               V_e = Mach * sqrt(gamma * T_e * R);
               rho_e = P_amb / (R * T_e);
           
            end 
       
       flowrate = Cd * rho_e * A_throat * V_e;
       
       % Compute thrust in phase 2
       Th = flowrate * V_e + (P_e - P_amb) * A_throat;
       % Compute drag 
       D = (rho_air / 2) .* (V_total).^2 * CD * A_bottle;
       
       % Velocity 
       %dx_dt -> change in the horizontal distance
       dx_dt = initial_data(3);
       % dz_dt -> change in the vertical distance
       dz_dt = initial_data(4);
       
       % Acceleration in x and z
       % dvx_dt -> change in the horizontal velocity with time
       dvx_dt = ((Th - D) * hx) ./ initial_data(7);
       % dvz_dt -> change in the vertical velocity with time
       dvz_dt = (((Th - D) * hz) - initial_data(7) * g)./ initial_data(7);
       
       % there is no change in the volume 
       dvol_dt = 0;
       % dmassa_dt -> change in the mass of air with time
       dmassa_dt = flowrate;
       % dmassr_dt -> change in the mass of the rocket with time
       dmassr_dt = -flowrate;
       
       % to return the derivatives as a column vector for the ode call
       dx = [dx_dt;dz_dt;dvx_dt;dvz_dt;dvol_dt;dmassa_dt;dmassr_dt];
       t2 = [t2 ; t];
   % --------------------------- Phase 3 --------------------------- %
   
       else
       
           V_total = sqrt((initial_data(3).^2)+(initial_data(4).^2));
           hx = initial_data(3)/V_total;
           hz = initial_data(4)/V_total;

           % There is no thrust on phase 3
           Th = 0;
           % Compute drag 
           D = (rho_air / 2) .* (V_total).^2 * CD * A_bottle;

           % Velocity 
           % dx_dt -> change in the horizontal distance
           dx_dt = initial_data(3);
           % dz_dt -> change in the vertical distance
           dz_dt = initial_data(4);

           % Acceleration in x and z
           % dvx_dt -> change in the horizontal velocity with time
           dvx_dt = ((Th - D) * hx) ./ initial_data(7);
           % dvz_dt -> change in the vertical velocity with time
           dvz_dt = (((Th - D) * hz) - initial_data(7) * g)./ initial_data(7);

           % there is no change in the volume 
           dvol_dt = 0;
           % there is no change in the mass of air with time
           dmassa_dt = 0;
           % there is no change in the mass of the rocket with time
           dmassr_dt = 0;

           dx = [dx_dt;dz_dt;dvx_dt;dvz_dt;dvol_dt;dmassa_dt;dmassr_dt];
           t3 = [t3 ; t];
           
       end
       
    end   
   
end 