function [air,plant_nodes,e_use] = unitary_system(building,air,T_supply,T_air,fan_avail,e_use,plant_nodes,i)
name = building.unitary_sys.fan_name{i};
f = f_index(name,building.fans.name);
air_flow = air.flow;
[air,P_fan] = fan_calc(building.fans,name,air,fan_avail(f));

switch building.unitary_sys.type{i}
    case 'FourPipeFanCoil'
        [air,plant_nodes,e_use] = unit_cool_coil(building,air,T_supply,T_air,e_use,plant_nodes,i);
        [air,plant_nodes,e_use] = unit_heat_coil(building,plant_nodes,air,T_supply,e_use,i);
    case 'PackagedTerminalAirConditioner'
        %% need to infer part-load-ratio (PLR), run on-cycle & off-cycle and sum power by PLR factor
        [air,plant_nodes,e_use] = unit_cool_coil(building,air,T_supply,T_air,e_use,plant_nodes,i);
        [air,plant_nodes,e_use] = unit_heat_coil(building,plant_nodes,air,T_supply,e_use,i);
    case 'UnitHeater'
        [air,plant_nodes,e_use] = unit_heat_coil(building,plant_nodes,air,T_supply,e_use,i);
end
if isfield(building.unitary_sys,'reheat_name') && ~isnan(building.unitary_sys.reheat_name(i)) && ~isempty(building.unitary_sys.reheat_name{i})
    disp('need to add dehumid and re-heat capability to unitary systems')
end
e_use.fan_elec = e_use.fan_elec + P_fan;
air.flow = air_flow;
end%Ends function unitary_system

function [air,plant_nodes,e_use] = unit_cool_coil(building,air,T_supply,T_air,e_use,plant_nodes,i)
name = building.unitary_sys.cool_coil{i};
switch building.unitary_sys.cool_coil_type{i}
    case {'Coil:Cooling:DX:TwoSpeed';'Coil:Cooling:DX:SingleSpeed'}
        [air,P_coil] = cooling_coil(building,name,air,T_air,T_supply);
        e_use.cool_elec = e_use.cool_elec + P_coil;
    case {'Coil:Cooling:Water:DetailedGeometry';'Coil:Cooling:Water'}
        [air,m_w,actual_load,out] = water_coil(building,name,air,T_supply,plant_nodes);
        plant_nodes.demand_flow(out) = m_w;%kg/s
        plant_nodes.load(out) = actual_load;%load in Watts (negative because air leaves hotter)
end
end%Ends function unit_cool_coil

function [air,plant_nodes,e_use] = unit_heat_coil(building,plant_nodes,air,T_supply,e_use,i)
name = building.unitary_sys.heat_coil{i};
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
end%Ends function unit_heat_coil
