function [air,plant_nodes,e_use] = unitary_heat_cool(building,fan_avail,air,T_supply,T_air,e_use,plant_nodes,i)
name = building.unitary_heat_cool.cool_coil{i};
if air.m_dot>1e-8
    if isfield(building.unitary_heat_cool,'cool_coil_type')
        switch building.unitary_heat_cool.cool_coil_type{i}
            case {'Coil:Cooling:DX:TwoSpeed';'Coil:Cooling:DX:SingleSpeed'}
                [air,P_coil] = cooling_coil(building,name,air,T_air,T_supply);
                e_use.cool_elec = e_use.cool_elec + P_coil;
            case {'Coil:Cooling:Water:DetailedGeometry';'Coil:Cooling:Water'}
                [air,m_w,actual_load,out] = water_coil(building,name,air,T_supply,plant_nodes);
                plant_nodes.demand_flow(out) = m_w;%kg/s
                plant_nodes.load(out) = actual_load;%load in Watts (negative because air leaves hotter)
        end
    end
    name = building.unitary_heat_cool.heat_coil{i};
    c = f_index(name,building.coils.heating.name);
    switch building.coils.heating.type{c}
        case {'Fuel';'Electric'}
            cp_supply = 1006 + 1860*air.w;
            load_heat = max(0,air.m_dot*cp_supply*(T_supply-air.T));
            if load_heat>0
                air.T = air.T + load_heat/(cp_supply*air.m_dot);
                air.h = psychometric(air,'h');
                if strcmp(building.coils.heating.type{c},'Fuel')
                    e_use.heat_gas = e_use.heat_gas + load_heat/building.coils.heating.efficiency(c); %heat added in W
                else
                    e_use.heat_elec = e_use.heat_elec + load_heat/building.coils.heating.efficiency(c); %heat added in W
                end
            end
        case {'Water'}
            [air,m_w,actual_load,out] = water_coil(building,name,air,T_supply,plant_nodes);
            plant_nodes.demand_flow(out) = m_w;%kg/s
            plant_nodes.load(out) = actual_load;%load in Watts (negative because air leaves hotter)
    end
end
name = building.unitary_heat_cool.fan_name{i};
f = f_index(name,building.fans.name);
[air,P_fan] = fan_calc(building.fans,name,air,fan_avail(f));
e_use.fan_elec = e_use.fan_elec + P_fan;
end%Ends function unitary_heat_cool