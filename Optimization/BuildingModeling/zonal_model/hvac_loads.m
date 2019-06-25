function [T,w,central_air,direct_air,mix,plenum,air_nodes,plant_nodes,frost,e_use,gains,window_gain,infiltration] = hvac_loads(building,schedules,air_nodes,T,w,plant_nodes,e_use,frost,weather,date,dt,des_day)
%Energy balance for multi-zone building (single step forward in time)

%%%Load values specific to this moment in time
[loads,gains,occupancy,mixing,infiltration,T_set,w_set,vect_2_sun,frost,T_mains,water_heat] = ...
    zone_loads(building,date,schedules,weather,T.zone',w.zone',dt,frost);
%%%Determine HVAC setpoints
[air_nodes,central_flow,direct_flow,mixed_air_loop,mixed_air_zone,mix,plenum,schedules,T,w,fan_avail,window_gain] = ...
    zone_ideal_HVAC(building,schedules,air_nodes,weather,mixing,infiltration,gains,T,w,T_set,w_set,occupancy',vect_2_sun,dt,des_day);
if isempty(des_day)
    %%%Compute performance of HVAC system 
    [central_air,e_use,plant_nodes] = sim_air_loops(building,air_nodes,mixed_air_loop,central_flow,w.loop,weather.DrybulbC,w.air,plant_nodes,e_use,fan_avail);
    [central_air,e_use,plant_nodes] = air_terminals(building,air_nodes,central_air,e_use,plant_nodes,T.central);
    [direct_air,plant_nodes,e_use] = unitary_systems(building,air_nodes,direct_flow,mixed_air_zone,T.direct,weather.DrybulbC,fan_avail,plant_nodes,e_use);%zone equipment not part of central air loops
    [plant_nodes,e_use,T.tank,T.tank_on] = manage_plant_equip(building,schedules,plant_nodes,e_use,dt,weather.DrybulbC,w.air,T_mains,T.tank,T.tank_on);
    e_use.water_gas = e_use.water_gas + water_heat/building.impact_factor.steam_efficiency;%unconnected water heat is considered purchased heat
end
end%Ends function hvac_loads

function [zone_air,e_use,plant_nodes] = air_terminals(building,air_nodes,zone_air,e_use,plant_nodes,terminal_T_supply)
%% re-heat or zone box cooling 
air_density = 1.204; %kg/m^3 %standard air density (1.204 kg/m3) adjusted for the local barometric pressure (standard barometric pressure corrected for altitude, ASHRAE 1997 HOF pg. 6.1).
for j = 1:1:length(building.hvac.terminals.name)
    switch building.hvac.terminals.type{j}
        case 'Uncontrolled'
            %%%
        case 'ConstantVolume:Reheat'
            %%%
        case 'VAV:Reheat'
            z = building.hvac.terminals.zone(j);
            m = building.zones.multiplier(z);
            if terminal_T_supply(z)>zone_air.T(z)
                air.T = zone_air.T(z);  air.w = zone_air.w(z);  air.m_dot = zone_air.m_dot(z)*m;  air.h = zone_air.h(z);%update intermediate air as it passes between equipment
                name = building.hvac.terminals.reheat_coil_name{j};
                c = f_index(name,building.coils.heating.name);
                switch building.coils.heating.type{c}
                    case {'Fuel';'Electric'}
                        cp_air = (1006 + 1860*air.w)*air_density;%Specific heat in J/m^3
                        Q_terminal = air.m_dot*cp_air*(terminal_T_supply(z) - air.T);
                        air.T = terminal_T_supply(z);
                        if strcmp(building.coils.heating.type{c},'Fuel')
                            e_use.heat_gas = e_use.heat_gas + Q_terminal/building.coils.heating.efficiency(c); %heat added in W
                        else
                            e_use.heat_elec = e_use.heat_elec + Q_terminal/building.coils.heating.efficiency(c); %heat added in W
                        end
                    case {'Water'}
                        [air,m_w,actual_load,out] = water_coil(building,name,air,terminal_T_supply(z),plant_nodes);
                        plant_nodes.demand_flow(out) = m_w;%kg/s
                        plant_nodes.load(out) = actual_load;%load in Watts (negative because air leaves hotter)
                end
                zone_air.T(z) = air.T; zone_air.w(z) = air.w; zone_air.m_dot(z) = air.m_dot/m; zone_air.h(z) = air.h;
            end
    end
end
end%Ends function air_terminals

function [direct_air,plant_nodes,e_use] = unitary_systems(building,air_nodes,direct_flow,mixed_air,T_direct,T_air,fan_avail,plant_nodes,e_use)
%% Unitary systems (zones not on an air loop)
direct_air = mixed_air;
for i = 1:1:length(building.unitary_sys.name)
    z = building.unitary_sys.zone(i);
    m = building.zones.multiplier(z);
    if z>0 && direct_flow(z)>0
        air.T = mixed_air.T(z,1);  air.w = mixed_air.w(z,1);  air.m_dot = mixed_air.m_dot(z,1)*m; air.flow = mixed_air.flow(z,1)*m; air.h = mixed_air.h(z,1);
        switch building.unitary_sys.type{i}
            case {'PackagedTerminalAirConditioner';'FourPipeFanCoil';'UnitHeater'}
                [air,plant_nodes,e_use] = unitary_system(building,air,T_direct(z),T_air,fan_avail,e_use,plant_nodes,i);
            case 'WindowAirConditioner'
                %not sure what to do about this (in zone, but not on air loop?)
            otherwise
                disp('need other types of unitary system')
        end
        direct_air.T(z) = air.T; direct_air.w(z) = air.w; direct_air.m_dot(z) = air.m_dot/m; direct_air.flow = air.flow/m; direct_air.h(z) = air.h;
    end
end
end%Ends function unitary_systems