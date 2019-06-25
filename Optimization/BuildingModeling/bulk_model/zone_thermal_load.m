function NetGain = zone_thermal_load(Build,Twall,Tdb,Teq,InternalGains,ExternalGains,dt)
%% net energy gain into the zone over the next time interval in kJ, assuming the zone goes to Teq

%estimate the temperature of the wall 1 step from now, T*_W in manual
f = .075;%fraction of wall R-value attributed to between wall and zone, remainder of R-value is between ambient and wall (needs to be 0.1 for small office)
Tsky = (16+273); %temperature in K that building is exhanging radiative heat transfer with
sigA = Build.VariableStruct.WallEmissivity*(Build.VariableStruct.WallArea*(1-Build.VariableStruct.WindowWallRatio) + Build.VariableStruct.RoofArea)*5.670367e-8/1000;% kW/K^4
UA_window = Build.VariableStruct.WallArea*Build.VariableStruct.WindowWallRatio*Build.VariableStruct.WindowUvalue/1000;
UA_wall = (Build.VariableStruct.WallArea*(1-Build.VariableStruct.WindowWallRatio)/Build.VariableStruct.WallRvalue + Build.VariableStruct.RoofArea/Build.VariableStruct.RoofRvalue)/1000;% Transmittance in kW/K, see the following about whole wall R-Value http://web.ornl.gov/sci/buildings/docs/Thermal-Performance-and-Wall-Ratings.pdf

Twall_est = Twall + ((Tdb - Twall)*(1/(1-f))*UA_wall + (Teq - Twall)*(1/f)*UA_wall + ExternalGains + (Tsky^4 - (273+Twall)^4)*sigA).*dt/(Build.VariableStruct.Capacitance*Build.Area);%estimate of wall temperature at t+1
Twall_est = 0.5*Twall + 0.25*Twall_est + 0.25*(Twall + ((Tdb - Twall_est)*(1/(1-f))*UA_wall + (Teq - Twall_est)*(1/f)*UA_wall + ExternalGains + (Tsky^4 - (273+Twall_est).^4)*sigA).*dt/(Build.VariableStruct.Capacitance*Build.Area));%average wall temp now with wall temp at t+1
NetGain = ((Tdb - Teq)*UA_window + (Twall_est - Teq)*(1/f)*UA_wall + InternalGains)*dt; 
end%Ends function zone_thermal_load