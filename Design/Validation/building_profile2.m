function [Cooling, Heating, FanPower,Damper,T_profile] = building_profile2(Build,Date,InternalGains,ExternalGains,Tdb,RH,T_init)
% function [Cooling, Heating, FanPower,Tzone,Twall,Damper] = building_profile(Build,Date,InternalGains,ExternalGains,Tdb,RH,Tzone_init,Twall_init)
% This function estimates the energy profile of a building (electric, heating and cooling kW)
% Build is a structure of building parameters
% Weather is an hourly weather profile (dry bulb, wet bulb, and relative humidity)
% Date is a vector of points in datenum format at which you are requesting the electric cooling & heating power

nS = length(Date); % number timesteps
m = length(T_init);
T_profile = zeros(nS+1,m);
Tsupply = zeros(nS,1);
T_profile(1,:) = T_init;
TsetC = load_sched(Build,[],Date,{'TsetC'},[]);
TsetH = load_sched(Build,[],Date,{'TsetH'},[]);
if nS == 1
    dt = 1;
else
    dt1 = Date(2) - Date(1);
    dt = (24*3600) * (Date - [Date(1)-dt1;Date(1:end-1)]); % duration of each time segment [seconds]
end
Cooling = zeros(nS,1);
Heating = zeros(nS,1);
AirFlow = zeros(nS,1);
Damper = zeros(nS,1);

R_total = 1/((Build.VariableStruct.WallArea*(1-Build.VariableStruct.WindowWallRatio)/Build.VariableStruct.WallRvalue + Build.VariableStruct.RoofArea/Build.VariableStruct.RoofRvalue)/1000);% Transmittance in kW/K, see the following about whole wall R-Value http://web.ornl.gov/sci/buildings/docs/Thermal-Performance-and-Wall-Ratings.pdf
R = [.3,.25,.25,.2]*R_total;
C = [10 [.2,.6,.2]*Build.VariableStruct.Capacitance]*Build.Area;
%% Heating/cooling dynamics
%need to step through time to account for moments with no heating/cooling
%as the building moves between the heat setpoint and the cooling setpoint
%note the conversion factor of 5.68 hft^2*F/Btu per 1 m^2*K/W
for t = 1:1:nS  
    if strcmp(Build.Name,'SmallOffice') && ~Build.VariableStruct.swOffPeakHvac && (TsetC(t) > 26 || TsetH(t) < 19) && (TsetC(max(1,t-round(3600/dt(t)))) > 26 || TsetH(max(1,t-round(3600/dt(t)))) < 19)
        %%HVAC is 'off'
        Tsupply(t) = nan;
        if T_profile(t,1)<21.5
            Heating(t) = 0.25*3.7382; %constant heating with no air flow (why does Eplus do this?)
        end
    else
        [Heating(t),Cooling(t),Damper(t),AirFlow(t),Tsupply(t)] = hvac_logic2(Build,R,C,T_profile(t,:),Tdb(t),RH(t),TsetC(t),TsetH(t),InternalGains(t),ExternalGains(t),dt(t));
    end
    [Heating(t),Cooling(t),T_profile(t+1,:)] = building_response2(Build,R,C,Tdb(t),RH(t),dt(t),InternalGains(t),ExternalGains(t),Cooling(t),Heating(t),Tsupply(t),AirFlow(t),Damper(t),T_profile(t,:));
%     [Heating(t),Cooling(t),Tzone(t+1),Twall(t+1)] = building_response(Build,Tdb(t),RH(t),dt(t),InternalGains(t),ExternalGains(t),Cooling(t),Heating(t),Tsupply(t),AirFlow(t),Damper(t),Tzone(t),Twall(t));
end
Cooling(abs(Cooling) < 1e-2) = 0;
FanPower = AirFlow*Build.VariableStruct.FanPower;
end%Ends function building_profile

function [Heating,Cooling,Damper,AirFlow,Tsupply] = hvac_logic2(Build,R,C,T_profile,Tdb,RH,TsetC,TsetH,InternalGains,ExternalGains,dt)
rho_Air = 1.225; % kg/m^3
% Fresh = true; %require fresh air recirculation
Damper = 1;
Heating = 0;
Cooling = 0;
Tsupply = nan;
%% Ambient dewpoint
P = 101.325; % atmospheric pressure (kPa)
Tdb_K = Tdb+273.15; %Tdb (Kelvin)
satP = exp((-5.8002206e3)./Tdb_K + 1.3914993 - 4.8640239e-2*Tdb_K + 4.1764768e-5*Tdb_K.^2 - 1.4452093e-8*Tdb_K.^3 + 6.5459673*log(Tdb_K))/1000; %saturated water vapor pressure ASHRAE 2013 fundamentals eq. 6 in kPa valid for 0 to 200C
P_H2O = RH/100.*satP; % kPa
Tdp = 6.54 + 14.526*log(P_H2O) + 0.7389*log(P_H2O).^2 + 0.09486*log(P_H2O).^3 + 0.4569*(P_H2O).^0.1984; %Dew point from partial pressure of water using ASHRAE 2013 Fundamentals eqn 39 valid from 0C to 93C
Cp_amb = 1.006 + 1.86*(.621945*(P_H2O./(P-P_H2O))); % kJ/kg*K
m_v_air = .621945*(P_H2O./(P-P_H2O));%mass fraction of water in air 

%% Specific heat of air
Tdp_set = 273+Build.VariableStruct.DPset;
P_H2O_dp = exp((-5.8002206e3)./Tdp_set + 1.3914993 - 4.8640239e-2*Tdp_set + 4.1764768e-5*Tdp_set.^2 - 1.4452093e-8*Tdp_set.^3 + 6.5459673*log(Tdp_set))/1000; %saturated water vapor pressure at dewpoint
m_v_set = .621945*(P_H2O_dp/(P-P_H2O_dp));%mass fraction of water in air at desired RH
Cp_build = 1.006 + 1.86*m_v_set; %kJ/kg
Latent_H2O = 2500.8 - 2.36*Build.VariableStruct.DPset +.0016*Build.VariableStruct.DPset^2 - 0.00006*Build.VariableStruct.DPset^3;% kJ/kg of latent heat
MinFlow = Build.VariableStruct.Volume*rho_Air*Build.VariableStruct.AirChangePerHr*(1/3600); %kg/s of flow
MaxFlow = 7.6*MinFlow;
Flow = MinFlow;
NetGain = zone_thermal_load2(Build,[TsetH,T_profile(2:end),Tdb],R,C,InternalGains,ExternalGains,dt);
if Build.VariableStruct.swFixFan %small office
    Flow = 10*MinFlow;
    Damper = .1;
end
if (min(TsetH,Build.VariableStruct.ColdAirSet) - Tdb)*Cp_amb*MinFlow*dt - NetGain > 0 
    %% Actively heating  
    if Build.VariableStruct.swFixFan %small office
        Tsupply =  TsetH - NetGain/(Flow*Cp_amb*dt); % temperature of supply air to provide this cooling
    elseif strcmp(Build.Name,'LO_DC_240_1zone_noFC') && TsetH < 19
        Damper = 0;
        Flow = MinFlow;
    else
        Tsupply = TsetH - NetGain/(Flow*Cp_amb*dt); % temperature of supply air to provide this heating
        Tsupply = min(65,Tsupply); % temperature of supply air to provide this heating 
        Flow = max(MinFlow,NetGain/((TsetH - Tsupply)*Cp_build*dt)); % mass flow of air to provide this heating
        Damper = MinFlow/Flow;
    end
else
    NetGain = zone_thermal_load2(Build,[TsetC,T_profile(2:end),Tdb],R,C,InternalGains,ExternalGains,dt);
    if (Tdb - TsetC)*Cp_amb*MaxFlow*dt + NetGain < 0
        %% Passive ventilation with exterior air
        if ~Build.VariableStruct.swFixFan
            if NetGain < 0 && TsetH < 19
                Damper = 0;
                Flow = MinFlow;
            elseif (Tdb - TsetC)*Cp_amb*MinFlow*dt + NetGain < 0
                Damper = min(1,(Build.VariableStruct.ColdAirSet - T_profile(1))/(Tdb - T_profile(1)));%if outdoor air is too cold, dilute with recirculation air
                Flow = max(MinFlow,.33*MinFlow/Damper);
            else
                Damper = 1;
                Flow = -NetGain/((Tdb - TsetC)*Cp_amb*dt);
            end
        end
    else
        %% Actively cooling
        if Build.VariableStruct.swFixFan %small office
            Tsupply =  TsetC - NetGain/(Flow*Cp_amb*dt); % temperature of supply air to provide this cooling
        else
            Tsupply = Build.VariableStruct.ColdAirSet;
            Flow = max(MinFlow, NetGain/((TsetC-Tsupply)*Cp_build*dt)); % mass flow of air to provide this cooling
            if Tdb < T_profile(1) % find economizer position
                Damper = max(MinFlow/Flow,min(1,(Tsupply - TsetC)./(Tdb - TsetC)));
            else
                Damper = MinFlow/Flow; % treat as little ouside air as possible
            end
        end 
    end
end
%% De-Humidification
if Build.VariableStruct.swDehumid && m_v_air>m_v_set %estimate de-humidification with Cp (need to add latent heat)
%             Cooling  = Cooling + Flow*Cp_Air*(Tmix - Build.VariableStruct.DPset)/n; %dehumidification energy in kJ/s
%             Heating  = Heating + Flow*Cp_Air*(Tmix - Build.VariableStruct.DPset)/n; %dehumidification energy in kJ/s
    Latent = Latent_H2O*Flow*Damper*max(0,m_v_air - m_v_set); %latent heat of condensation in kJ/s
    Cooling  = Cooling + Latent;
end
AirFlow = Flow/rho_Air; % flow rate in m^3/s
end%Ends function hvac_logic

function [Heating,Cooling,T_final] = building_response2(Build,R,C,Tdb,RH,dt,InternalGains,ExternalGains,Cooling,Heating,Tsupply,AirFlow,Damper,T_init)
%%Note this function is very similar to, but cannot be replaced by building_simulate
%%This finds the heating to achieve a Tsupply, where building_simulate
%%already knows the Heating and Cooling
%% there are m + 1 states where the wall is divided into m segments (+1 zone air temperature)
rho_Air = 1.225; % kg/m^3
n = floor(dt/30);% 30 second steps size for smotth simulation of temperature profile
%% Ambient dewpoint
P = 101.325; % atmospheric pressure (kPa)
Tdb_K = Tdb+273.15; %Tdb (Kelvin)
satP = exp((-5.8002206e3)./Tdb_K + 1.3914993 - 4.8640239e-2*Tdb_K + 4.1764768e-5*Tdb_K.^2 - 1.4452093e-8*Tdb_K.^3 + 6.5459673*log(Tdb_K))/1000; %saturated water vapor pressure ASHRAE 2013 fundamentals eq. 6 in kPa valid for 0 to 200C
P_H2O = RH/100.*satP; % kPa
Cp_amb = 1.006 + 1.86*(.621945*(P_H2O./(P-P_H2O))); % kJ/kg*K
%% Specific heat of air
Tdp_set = 273+Build.VariableStruct.DPset;
P_H2O_dp = exp((-5.8002206e3)./Tdp_set + 1.3914993 - 4.8640239e-2*Tdp_set + 4.1764768e-5*Tdp_set.^2 - 1.4452093e-8*Tdp_set.^3 + 6.5459673*log(Tdp_set))/1000; %saturated water vapor pressure at dewpoint
m_v_set = .621945*(P_H2O_dp/(P-P_H2O_dp));%mass fraction of water in air at desired RH
Cp_build = 1.006 + 1.86*m_v_set; %kJ/kg
%%
Tsky = (16+273); %temperature in K that building is exhanging radiative heat transfer with
sigA = Build.VariableStruct.WallEmissivity*(Build.VariableStruct.WallArea*(1-Build.VariableStruct.WindowWallRatio) + Build.VariableStruct.RoofArea)*5.670367e-8/1000;% kW/K^4
UA_window = Build.VariableStruct.WallArea*Build.VariableStruct.WindowWallRatio*Build.VariableStruct.WindowUvalue/1000;
m = length(T_init);
T = zeros(n+1,m+1);
T(1,1:m) = T_init;
T(:,end) = Tdb;

Flow = AirFlow*rho_Air; % flow rate in m^3/s
if isnan(Tsupply)
    mode='passive';
else
    mode = 'active';
end
for i = 1:1:n
    Tmix = Damper*T(i,m+1) + (1-Damper)*T(i,1);%note T(i,m+1) is Tdb and T(i,1) = T_air
    Cp_Air = Damper*Cp_amb + (1-Damper)*Cp_build;
    if strcmp(mode,'passive')
        Tsupply = Tmix;
    else
        Q = Cp_Air*(Tsupply-Tmix)*Flow/n;
        if Q>0
            Heating = Heating + Q;
        else
            Cooling = Cooling - Q;
        end
    end
    if Flow ==0
        T(i+1,1) = T(i,1) + ((T(i,m+1) - T(i,1))*UA_window + (T(i,2) - T(i,1))/R(1) + InternalGains + Heating)*dt/(n*C(1)); %net change in temperature
    else
        T(i+1,1) = T(i,1) + ((T(i,m+1) - T(i,1))*UA_window + (T(i,2) - T(i,1))/R(1) + InternalGains + (Tsupply - T(i,1))*Cp_Air*Flow)*dt/(n*C(1)); %net change in temperature
    end
    for j = 2:m-1
        T(i+1,j) = T(i,j) + ((T(i,j+1) - T(i,j))/R(j) + (T(i+1,j-1) - T(i,j))/R(j-1)).*dt/(n*C(j));
    end
    T(i+1,m) = T(i,m) + ((T(i,m+1) - T(i,m))/R(m) + (T(i+1,m-1) - T(i,m))/R(m-1) + ExternalGains + (Tsky^4 - (273+T(i,m))^4)*sigA).*dt/(n*C(m));
end
T_final = T(end,1:m);
end%Ends function building_response2

function NetGain = zone_thermal_load2(Build,T,R,C,InternalGains,ExternalGains,dt)
%estimate the temperature of the wall 1 atep from now, T*_W in manual
Tsky = (16+273); %temperature in K that building is exhanging radiative heat transfer with
sigA = Build.VariableStruct.WallEmissivity*(Build.VariableStruct.WallArea*(1-Build.VariableStruct.WindowWallRatio) + Build.VariableStruct.RoofArea)*5.670367e-8/1000;% kW/K^4
UA_window = Build.VariableStruct.WallArea*Build.VariableStruct.WindowWallRatio*Build.VariableStruct.WindowUvalue/1000;
m = length(T)-1;

T_est = zeros(1,m);
for j = 2:m-1
    T_est(j) = T(j) + ((T(j+1) - T(j))/R(j) + (T(j-1) - T(j))/R(j-1)).*dt/C(j);
end
T_est(m) = T(m) + ((T(m+1) - T(m))/R(m) + (T_est(m-1) - T(m))/R(m-1) + ExternalGains + (Tsky^4 - (273+T(m))^4)*sigA).*dt/C(m);
for j = 2:1:m-1
T_est(j) = 0.5*T(j) + 0.25*T_est(j) + 0.25*(T(j) + ((T(j+1) - T_est(j))/R(j) + (T(j-1) - T_est(j))/R(j-1)).*dt/C(j));
end
T_est(m) = 0.5*T(m) + 0.25*T_est(m) + 0.25*(T(m) + ((T(m+1) - T_est(m))/R(m) + (T(m-1) - T_est(m))/R(m-1) + ExternalGains + (Tsky^4 - (273+T_est(m)).^4)*sigA).*dt/C(m));%average wall temp now with wall temp at t+1

NetGain = ((T(end) - T(1))*UA_window + (T_est(2) - T(1))/R(1) + InternalGains)*dt; %net energy gain into the zone in kJ, assuming the zone goes to TsetH
end%Ends function zone_thermal_load2