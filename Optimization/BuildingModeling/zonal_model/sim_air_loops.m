function [zone_air,e_use,plant_nodes] = sim_air_loops(building,air_nodes,supply_air,central_flow,w_loop,T_air,m_v_air,plant_nodes,e_use,fan_avail)
for l = 1:1:length(supply_air.T)  %loop through equipment on loop
    on_loop = f_index(l,building.air_supply_equip.loop);
    air.T = supply_air.T(l);air.w = supply_air.w(l);air.m_dot = supply_air.m_dot(l);air.h = supply_air.h(l);%update intermediate air as it passes between equipment
    for j = 1:1:length(on_loop)
        %% making sure simulation of these components is in the order they are on the supply branch
        e = f_index(j,building.air_supply_equip.branch_order(on_loop));
        k = on_loop(e);
        type = building.air_supply_equip.type{k};
        name = building.air_supply_equip.name{k};
        out = building.air_supply_equip.outlet{k};
        switch type
            case 'AirLoopHVAC:OutdoorAirSystem'
                %do nothing, already mixed air
            case {'Coil:Cooling:Water:DetailedGeometry';'Coil:Cooling:Water';'Coil:Heating:Water'}
                [air,m_w,actual_load,out] = water_coil(building,name,air,air_nodes.supply_T_set(out),plant_nodes);
                plant_nodes.demand_flow(out) = m_w;%kg/s
                plant_nodes.load(out) = actual_load;%load in Watts (negative for heating coil because air leaves hotter)
                if strcmp(type,'Coil:Cooling:Water') && actual_load<0
                    disp('wtf')
                end
            case {'CoilSystem:Cooling:DX';'DX:TwoSpeed';'DX:SingleSpeed'}
                T_wb = psychometric('T',T_air,'w',m_v_air,'Twb');
                [air,P_coil] = cooling_coil(building,name,air,T_air,air_nodes.supply_T_set(out));
                e_use.cool_elec = e_use.cool_elec + P_coil;
            case {'Coil:Heating:Fuel';'Coil:Heating:Electric'}
                spec_heat_supply = 1006 + 1860*air.w; %J/kg*K
                Q_treat = air.m_dot*spec_heat_supply*(air_nodes.supply_T_set(out)-air.T);
                c = f_index(name,building.coils.heating.name);
                if Q_treat>1
                    air.T = air_nodes.supply_T_set(out);
                    air.h = air.h + Q_treat/air.m_dot;
                    if strcmp(type,'Coil:Heating:Fuel')
                        e_use.heat_gas = e_use.heat_gas + Q_treat/building.coils.heating.efficiency(c); %heat added in W
                    else
                        e_use.heat_elec = e_use.heat_elec + Q_treat/building.coils.heating.efficiency(c); %heat added in W
                    end
                end
            case {'Fan:ConstantVolume';'Fan:VariableVolume';'Fan:OnOff';'Fan:SystemModel'}
                f = f_index(name,building.fans.name);
                [air,P_fan] = fan_calc(building.fans,name,air,fan_avail(f));
                e_use.fan_elec = e_use.fan_elec + P_fan;
            case 'Humidifier:Steam:Electric'
                if w_loop(l)>air.w %see page 1250
                    %steam specific enthalpy = 2676125 J/kg
                    den_water = 995;%kg/m^3
                    m_w = (w_loop(l)-air.w)*air.m_dot;
                    air.h = air.h + m_w/air.m_dot*2676125;
                    air.T = psychometric(air,'T');
                    %% need to check it is not over saturated
                    k = f_index(name,building.humidifiers.name);
                    e_use.heat_elec = e_use.heat_elec + m_w/(den_water*building.humidifiers.flow_capacity(k))*building.humidifiers.max_power(k) + building.humidifiers.fan_power(k) + building.humidifiers.standby_power(k);
                end
            case 'AirLoopHVAC:UnitaryHeatCool'
                k = f_index(name,building.unitary_heat_cool.name);
                [air,plant_nodes,e_use] = unitary_heat_cool(building,fan_avail,air,air_nodes.supply_T_set(out),T_air,e_use,plant_nodes,k);
            otherwise
                disp('need to add components in sim_air_loops');
        end 
    end
    supply_air.T(l) = air.T; supply_air.w(l) = air.w; supply_air.m_dot(l) = air.m_dot; supply_air.h(l) = air.h;%update intermediate air as it passes between equipment
end

air_density = 1.204; %kg/m^3 %standard air density (1.204 kg/m3) adjusted for the local barometric pressure (standard barometric pressure corrected for altitude, ASHRAE 1997 HOF pg. 6.1).
n = length(building.zones.name);
zone_air.T = zeros(n,1);
zone_air.w = zeros(n,1);
zone_air.h = zeros(n,1);
zone_air.flow = central_flow;
zone_air.m_dot = air_density*central_flow;
for l = 1:1:length(supply_air.T) 
    A = (building.hvac.zone_2_loop(l,:)>0);
    zone_air.T(A) = supply_air.T(l);
    zone_air.w(A) = supply_air.w(l);
    zone_air.h(A) = supply_air.h(l);
end

%% exhaust fan power
for i = 1:1:length(building.fans.name)
    if building.fans.exhaust(i) && ~isnan(building.fans.flow_rate(i))
        fan_flow = building.fans.flow_rate(i)*fan_avail(i);%volumetric flow (m^3/s)  %See page 1601 of the reference manual
        e_use.fan_elec = e_use.fan_elec + fan_flow.*building.fans.pressure_rise(i)/(building.fans.fan_efficiency(i));% m3/s * Pa = N*m/s = W
    end
end
end%Ends function sim_air_loops