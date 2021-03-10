% Andrew Logue
% 12/6/2020
% Project 2 - ASEN 2012

function [dx] = fun(~,x_i)
    % function for calculating rocket trajectory
    % input: initial conditions vector x_i in the bellow format:
    % x_i(1)- distance (x position)
    % x_i(2)- height (z position)
    % x_i(3)- velocity (x component)
    % x_i(4)- velocity (z component)
    % x_i(5)- volume of air
    % x_i(6)- mass of air
    % x_i(7)- mass of rocket

    % output: derivative vector dx in the bellow format:
    % dx(1)- horizontal velocity
    % dx(2)- vertical velocity
    % dx(3)- horizontal acceleration
    % dx(4)- vertical acceleration
    % dx(5)- change of volume of air in bottle rocket
    % dx(6)- change of mass of air in bottle rocket
    % dx(7)- change of mass of the bottle rocket

    global emptyBottleVolume standL theta totalPressure R g Cd P_f...
        waterDensity airDensity CD bottleArea throatArea airMass_0...
        atmosphericPressure airVolume_0 specificHeatRatio z_0


    % PHASE 1
    if emptyBottleVolume > x_i(5)
         % check if liftoff occured -> rocket would no longer be on test stand
        % pythag
        if sqrt((x_i(1)^2) + (x_i(2) - z_0)^2) > standL
            % has launched off the test stand
            x_1 = x_i(3) / sqrt((x_i(4).^2) + (x_i(3).^2));
            z_1 = x_i(4) / sqrt((x_i(4).^2) + (x_i(3).^2));

        else
            % still on the test stand
            x_1 = cos(theta);
            z_1 = sin(theta);

        end
        
        % forces on rocket: pressure, thrust, drag
        pressure = totalPressure.*((airVolume_0./x_i(5)).^...
            specificHeatRatio);
        thrust = 2.*Cd.*(pressure - atmosphericPressure).*throatArea;
        drag = (CD * bottleArea * (sqrt((x_i(4).^2) + (x_i(3).^2)))).*...
            (airDensity/2);
        
        % calculate dx derivaties
        % dx(1)- horizontal velocity
        % dx(2)- vertical velocity
        % dx(3)- horizontal acceleration
        % dx(4)- vertical acceleration
        % dx(5)- change of volume of air in bottle rocket
        % dx(6)- change of mass of air in bottle rocket (none = 0)
        % dx(7)- change of mass of the bottle rocket
        dx = [x_i(3); x_i(4); ((thrust - drag) * x_1) ./ x_i(7);...
            (((thrust - drag) * z_1) - x_i(7) * g) ./ x_i(7);...
            throatArea * Cd * sqrt((2 / waterDensity) * (totalPressure *...
            ((airVolume_0 / x_i(5))^specificHeatRatio) -...
            atmosphericPressure)); 0; throatArea .* (-Cd) .*...
            sqrt(2.*waterDensity.*(pressure - atmosphericPressure))];

    elseif emptyBottleVolume <= x_i(5)   
        % has launched off the test stand
        x_1 = x_i(3) / sqrt((x_i(4).^2) + (x_i(3).^2));
        z_1 = x_i(4) / sqrt((x_i(4).^2) + (x_i(3).^2));
        pressure =  P_f * (x_i(6)/airMass_0)^specificHeatRatio;

            % Check condition that separates the 2nd phase from the 3rd
            % (preassure) 
            if atmosphericPressure < pressure
                % calculate critical pressure
                criticalPressure = (pressure)*(2./(specificHeatRatio +...
                    1)).^(specificHeatRatio / (specificHeatRatio-1));
                airTemp = pressure / ((x_i(6) / emptyBottleVolume)*R);

                % calculate final temp, density, velocity, and pressure
                if atmosphericPressure < criticalPressure
                    P_f2 = criticalPressure;
                    T_f = airTemp * (2/(specificHeatRatio+1));
                    v_f = sqrt(T_f*R*specificHeatRatio);
                    density_f = criticalPressure / (R*T_f);

                else
                    % calculate Mach
                    Mach = sqrt(((pressure/atmosphericPressure)^...
                        (((specificHeatRatio-1)/specificHeatRatio))-1)*...
                        (2/(specificHeatRatio-1)));
                    P_f2 = atmosphericPressure;
                    T_f = airTemp / (1+((specificHeatRatio-1)/2) * Mach^2);
                    v_f = Mach * sqrt(T_f*R*specificHeatRatio);
                    density_f = atmosphericPressure / (R * T_f);

                 end

            % Compute the air's mass flow rate
            airflow = Cd * density_f * throatArea * v_f;

            % forces on rocket: pressure, thrust, drag
            % pressure already calculated
            thrust = airflow*v_f + (P_f2 - atmosphericPressure)*throatArea;
            drag = (airDensity / 2).*(sqrt((x_i(4).^2)+(x_i(3).^2))).^2 ...
                * CD * bottleArea;

            % calculate dx derivaties
            % dx(1)- horizontal velocity
            % dx(2)- vertical velocity
            % dx(3)- horizontal acceleration
            % dx(4)- vertical acceleration
            % dx(5)- change of volume of air in bottle rocket (none = 0)
            % dx(6)- change of mass of air in bottle rocket 
            % dx(7)- change of mass of the bottle rocket 
            dx = [x_i(3); x_i(4); ((thrust - drag) * x_1) ./ x_i(7);...
                (((thrust - drag) * z_1) - (x_i(7) * g)) ./ x_i(7); 0;...
                -airflow; -airflow];
            
            else 
                % has launched off the test stand
                x_1 = x_i(3) / sqrt((x_i(4).^2) + (x_i(3).^2));
                z_1 = x_i(4) / sqrt((x_i(4).^2) + (x_i(3).^2));
                % PHASE 3
                % forces on rocket: pressure, thrust, drag
                % pressure already calculated
                % no thrust
                drag = (airDensity / 2).*(sqrt((x_i(4).^2)+(x_i(3).^2)))...
                    .^2 * CD * bottleArea;

                % calculate dx derivaties
                % dx(1)- horizontal velocity
                % dx(2)- vertical velocity
                % dx(3)- horizontal acceleration
                % dx(4)- vertical acceleration
                % dx(5)- change of volume of air in bottle rocket(none = 0)
                % dx(6)- change of mass of air in bottle rocket (none = 0)
                % dx(7)- change of mass of the bottle rocket (none = 0)
                dx = [x_i(3); x_i(4); (-drag * x_1) ./ x_i(7);...
                    ((-drag * z_1) - (x_i(7) * g)) ./ x_i(7); 0; 0; 0];
            end
    end

end
