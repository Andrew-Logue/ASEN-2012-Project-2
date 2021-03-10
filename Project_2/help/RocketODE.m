function [ derivatives ] = RocketODE(Time,States)
%
%
global TestStandLength Theta Pgage Pamb Cd ThroatArea CD BottleArea Rhoairamb...
    RhoWater Volbottle y0 VAirInit GammaGas g  MassAirInit R t2 t1 t3
%
%
%
% --------- (States In order)------------------
% 1- Mass of rocket;
% 2- Mass of Air
% 3- Volume of Air;
% 4- Velocity x;
% 5- Velocity z;
% 6- Range (X location);
% 7- Height (Z location);


% 

% PREDEFINE Constants:

%% Phase 1: 

if States(5) < Volbottle
    

    % Check if we still on stand or left:
    
    h = sqrt((States(1)^2)+(States(2)-y0)^2);
    
    if h <= TestStandLength
        TotalVeloc = sqrt( (States(4).^2) + (States(3).^2) );
        HeadingX = cosd(Theta);
        HeadingZ = sind(Theta);
    else
        TotalVeloc = sqrt( (States(4).^2) + (States(3).^2) );
        HeadingX = States(3)/TotalVeloc;
        HeadingZ = States(4)/TotalVeloc;
    end

    Pressure = ( ( VAirInit ./ States(5) ) .^ GammaGas ) .* (Pgage+Pamb) ; 
    Thrust = 2.* Cd .* ThroatArea .* ( Pressure - Pamb) ;
    Drag = ( Rhoairamb / 2) .* (TotalVeloc).^2 * CD*BottleArea; 

    % Define Derivatives

    % Accelration in x and Z
    dadt_X = ( (Thrust - Drag) * HeadingX) ./ States(7) ;
    dadt_Z =  ( ((Thrust - Drag) * HeadingZ) - States(7)*g ) ./ States(7) ;

    %Volume (how volume changes with time)
    DVolume_Dt = Cd * ThroatArea * sqrt ( (2/RhoWater) * ( ( (Pgage+Pamb) * (( VAirInit/States(5) ) ^ (GammaGas))) - Pamb ));

    %Mass (how mass of Rocket changes with time, it's )
    DMass_Dt = - Cd .* ThroatArea .* sqrt ( 2.*RhoWater.* ( Pressure - Pamb ) );

    derivatives = [ States(3) ; States(4) ; dadt_X; dadt_Z; DVolume_Dt ; 0 ; DMass_Dt] ;
%     derivatives = [ DMass_Dt; 0; DVolume_Dt; dadt_X; dadt_Z; States(3) ; States(4) ] ;
    t1 = [ t1 ; [Time Thrust] ];
%% Phase 2:

elseif States(5)>= Volbottle
%T and P of end states
% Tend = TAirInit * (( VAirInit/Volbottle) ^ (GammaGas-1) );
Pend = (Pgage+Pamb) * (( VAirInit/Volbottle) ^ (GammaGas) );

    TotalVeloc = sqrt( (States(4).^2) + (States(3).^2) );
    HeadingX = States(3)/TotalVeloc;
    HeadingZ = States(4)/TotalVeloc;
    
PressureCond = Pend * (States(6)/MassAirInit)^(GammaGas) ;

if PressureCond>Pamb
Density = States(6) / Volbottle;
Temp = PressureCond/(Density*R);
CriticalP = (PressureCond) * (2./(GammaGas+1)).^(GammaGas/(GammaGas-1));
  
if CriticalP > Pamb
    
%     Mach  = 1;
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
MassRocketFlowRate = -MassAirFlowRate;

Thrust = MassAirFlowRate *Vexit + (Pexit-Pamb)*ThroatArea ;
Drag = ( Rhoairamb / 2) .* (TotalVeloc).^2 * CD*BottleArea; 

dadt_X = ( (Thrust - Drag) * HeadingX) ./ States(7) ;
dadt_Z =  ( ((Thrust - Drag) * HeadingZ) - States(7)*g ) ./ States(7) ;

derivatives = [ States(3) ; States(4) ; dadt_X; dadt_Z; 0 ; -MassAirFlowRate ; MassRocketFlowRate] ;
% derivatives = [ MassRocketFlowRate; -MassAirFlowRate; 0; dadt_X; dadt_Z; States(3) ; States(4) ] ;
t2 = [ t2 ; Time ];


%% Phase 3: 

else

TotalVeloc = sqrt( (States(4).^2) + (States(3).^2) );
HeadingX = States(3)/TotalVeloc;
HeadingZ = States(4)/TotalVeloc;

Thrust = 0 ;
Drag = ( Rhoairamb / 2) .* (TotalVeloc).^2 * CD*BottleArea; 

dadt_X = ( (Thrust - Drag) * HeadingX) ./ States(7) ;
dadt_Z =  ( ((Thrust - Drag) * HeadingZ) - States(7)*g ) ./ States(7) ;

derivatives = [ States(3) ; States(4); dadt_X; dadt_Z; 0; 0 ; 0 ] ;
% derivatives = [ 0; 0; 0; dadt_X; dadt_Z; States(3) ; States(4) ] ;
t3 = [ t3 ; Time ];

end



end

