% Andrew Logue
% 12/6/2020
% Project 2 - ASEN 2012

function [thrust_1, thrust_2, thrust_3] = thrust(x)
    % function for calculating rocket trajectory
    % input: time vector t, and x vector from ODE45

    % output:
    % thrust_1- phase 1 thrust
    % thrust_2- phase 2 thrust
    % thrust_3- phase 3 thrust
    
    % include global variables
    global emptyBottleVolume throatArea airVolume_0 atmosphericPressure...
        gagePressure R Cd specificHeatRatio  airMass_0 P_f
    
    % get rows of x
    [r,~] = size(x);
    % initialize counter variables
    count_1 = 0;
    count_2 = 0;
    totalPressure = gagePressure + atmosphericPressure;
    % for each row
    
    
    for i = 1:r
        % if in phase 1
        if  emptyBottleVolume > x(i,5) 
            P(i) = ((airVolume_0 ./ x(i,5)).^specificHeatRatio).*(totalPressure);
            thrust_1(i) = 2.*Cd.*throatArea.*(P(i) - atmosphericPressure);

        elseif x(i,5) >= emptyBottleVolume
            % if in phase 2
            pressure = P_f * (x(i,6) / airMass_0)^specificHeatRatio;
            
            if atmosphericPressure < pressure
                
                D = x(i,6) / emptyBottleVolume;
                airTemp = pressure / (D * R);
                criticalPressure = (pressure) * (2./...
                    (specificHeatRatio+1)).^(specificHeatRatio /...
                    (specificHeatRatio - 1));
   
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
                    density_f = atmosphericPressure / (R*T_f);

                 end

                % Compute flow rate of air
                airflow = Cd * density_f * throatArea * v_f;
%                 count_1 ++
                count_1 = count_1+1;
                % get thrust in phase 2
                thrust_2(count_1) = airflow * v_f + (P_f2 - ...
                    atmosphericPressure) * throatArea;

            else
                % count_2 ++
                count_2 = count_2+1;
                % thrust in phase 3 = 0
                thrust_3(count_2) = 0;

            end
        end
    end
end

