function HVAC_out = manage_air_loops(building,profile,zone_flow,T_supply_zone,T_loop_supply,mixed_air,non_loop_mixed_air,Q_terminal,T_air)
reheat_zone = Q_terminal.*(Q_terminal>0); %re-heat or zone box cooling (need 0 for passive condition when Q_ideal_hvac == 0)
reheat_displaced = 0;%building.HVAC.zone_2_loop*(reheat_zone.*building.zones.multiplier);
[load_cool,load_heat] = system_treatment(building.HVAC,mixed_air,T_loop_supply,profile,reheat_displaced);
%need to update building.plant.loop.hot_water_supply_temperature & building.plant.loop.hot_water_supply_temperature with actual loop water temperature simulation so the plant water temperature is solved
[supply_air,e_use,fans] = sim_supply_equipment(building,mixed_air,T_air,profile,load_cool,load_heat);
[zone_air,e_use] = sim_zone_equipment(building,supply_air,zone_flow,reheat_zone,non_loop_mixed_air,T_air,T_supply_zone,e_use);%re-visit T_set_loop to make T_set_zone if reheat was subtracted from Q_loop

%% Add Plant loop to convert heating.water and cooling.water into electricity
e_use  = manage_plant_equip(building,e_use);

T_supply = zone_air.T';
m_v_supply = zone_air.w;
T_supply(zone_flow==0) =  T_supply_zone(zone_flow==0);%T_supply_loop(1);
m_v_supply(zone_flow==0) = mixed_air.w(1);
HVAC_out = {T_supply;m_v_supply;e_use.heat_elec;e_use.heat_fuel;e_use.cool_elec;fans};
end%Ends function manage_air_loops

function [load_cool,load_heat] = system_treatment(hvac,mixed_air,T_supply,profile,reheat_displaced)
%This is the "unitary system" part of the air loop
Q_fan = zeros(length(mixed_air.T),1);
for i = 1:1:length(hvac.loop.name)
    for j = 1:1:length(hvac.loop.component_name{i})
        switch hvac.loop.component_type{i}{j}
            case {'Fan:ConstantVolume';'Fan:VariableVolume';'Fan:OnOff'}
                air_est.T = T_supply(i);air_est.w = mixed_air.w(i);air_est.m_dot = mixed_air.m_dot(i);
                air_est.h = psychometric_h(air_est);
                f = nonzeros((1:length(hvac.fans.name))'.*strcmp(hvac.loop.component_name{i}{j},hvac.fans.name));
                [~,~,Q_to_air] = fan_calc(hvac.fans,hvac.loop.component_name{i}{j},'operation',air_est,profile.(hvac.fans.schedule{f}));
                Q_fan(i) = Q_fan(i) + Q_to_air;
        end
    end
end

spec_heat_supply = 1006 + 1860*mixed_air.w; %J/kg*K
%cooling in zone calculation is m*Cp*dT, not h2-h1, so need to ensure sensible load = m*Cp*(Tsupply -Tzone)
T_treated = T_supply -(reheat_displaced + Q_fan)./(mixed_air.m_dot.*spec_heat_supply);
A = (mixed_air.T>T_treated & (mixed_air.T-T_supply)<1e-2) | mixed_air.m_dot<1e-8; %was heating to get to T_set.heat, now would need cooling, consider as passive
T_treated(A) = mixed_air.T(A);
Q_treat = mixed_air.m_dot.*spec_heat_supply.*(T_treated-mixed_air.T);
%%Need to add de-humidifying case
load_cool = max(0,-Q_treat);
load_heat = max(0,Q_treat);
end%Ends function system_treatment

function [supply_air,e_use,fans] = sim_supply_equipment(building,mixed_air,T_air,profile,load_cool,load_heat)
%This is the "unitary system" part of the air loop
%loop through equipment on loop
hvac = building.HVAC;
loops = length(hvac.zone_2_loop(:,1));
e_use.cool_elec = 0;
e_use.heat_fuel = 0;
e_use.heat_elec = 0;
e_use.heat_water.return_flow = [];
e_use.heat_water.return_temp = [];
e_use.heat_water.plant_loop = [];
e_use.cool_water.return_flow = [];
e_use.cool_water.return_temp = [];
e_use.cool_water.plant_loop = [];

fans = zeros(1,length(hvac.fans.name));
for i = 1:1:loops
    if mixed_air.m_dot(i)<1e-8
        supply_air.T(i) = mixed_air.T(i);supply_air.w(i) = mixed_air.w(i);supply_air.m_dot(i) = mixed_air.m_dot(i);supply_air.h(i) = mixed_air.h(i);
    else
        intermediate_air.T = mixed_air.T(i);intermediate_air.w = mixed_air.w(i);intermediate_air.m_dot = mixed_air.m_dot(i);intermediate_air.h = mixed_air.h(i);%update intermediate air as it passes between equipment

        %%need to add controller so it leaves at T_target after heater & fans.
        %%heater load reduced by heat into fans. See bit about setpoint
        %%controllers (pg 1070 for heating, page 1064 for cooling)
        %% converge to a tolerance on sensible heating supply
        %first determine heating or cooling
        %if cooling, run full load cooling coil, then zero load cooling coil
        %and calculate the part-load-fraction, similar for heating
        for j = 1:1:length(hvac.loop.component_name{i})            
            type = hvac.loop.component_type{i}{j};
            name = hvac.loop.component_name{i}{j};
            switch type
                case 'AirLoopHVAC:OutdoorAirSystem'
                    %do nothing
                case 'AirLoopHVAC:UnitaryHeatCool'
                    k = nonzeros((1:length(building.HVAC.unitary_sys.name))'.*strcmp(name,building.HVAC.unitary_sys.name));
                    [intermediate_air,e_use] = unitary_sys_operate(intermediate_air,load_heat(i),load_cool(i),T_air,e_use,building.HVAC,k);
                case {'Coil:Cooling:Water:DetailedGeometry';'Coil:Cooling:Water'}
                    pl = building.plant.components.loop(nonzeros((1:length(building.plant.components.name))'.*strcmpi(name,building.plant.components.name)));%plant loop
                    [intermediate_air,water] = cooling_coil(hvac.coils.cooling,name,'operation',intermediate_air,load_cool(i),building.plant.loop.cold_water_supply_temperature(pl));
                    e_use.cool_water.plant_loop(end+1) = pl;
                    e_use.cool_water.return_flow(end+1) = water.m_dot;
                    e_use.cool_water.return_temp(end+1) = water.T;
                case {'CoilSystem:Cooling:DX';'DX:TwoSpeed';'DX:SingleSpeed'}
                    if ~any(strcmpi(name,hvac.coils.cooling.name))
                        name = hvac.coil_system.coil_name{nonzeros((1:length(hvac.coil_system.name))'.*strcmpi(name,hvac.coil_system.name))};
                    end
                    [intermediate_air,P_coil] = cooling_coil(hvac.coils.cooling,name,'operation',hvac.curves,intermediate_air,T_air,load_cool(i),intermediate_air.w);
                    e_use.cool_elec = e_use.cool_elec + P_coil;
                case 'Coil:Heating:Fuel'
                    [intermediate_air,Q_coil] = heating_coil(hvac.coils.heating,name,'operation',intermediate_air,load_heat(i));
                    e_use.heat_fuel = e_use.heat_fuel + Q_coil;
                case 'Coil:Heating:Electric'
                    [intermediate_air,Q_coil] = heating_coil(hvac.coils.heating,name,'operation',intermediate_air,load_heat(i));
                    e_use.heat_elec = e_use.heat_elec + Q_coil;
                case 'Coil:Heating:Water'
                    pl = building.plant.components.loop(nonzeros((1:length(building.plant.components.name))'.*strcmpi(name,building.plant.components.name)));%plant loop
                    [intermediate_air,water] = heating_coil(hvac.coils.heating,name,'operation',intermediate_air,load_heat(i),building.plant.loop.hot_water_supply_temperature(pl));
                    e_use.heat_water.plant_loop(end+1) = pl;
                    e_use.heat_water.return_flow(end+1) = water.m_dot;
                    e_use.heat_water.return_temp(end+1) = water.T;
                case {'Fan:ConstantVolume';'Fan:VariableVolume';'Fan:OnOff'}
                    f = nonzeros((1:length(hvac.fans.name))'.*strcmp(name,hvac.fans.name));
                    [intermediate_air,fans(f)] = fan_calc(hvac.fans,name,'operation',intermediate_air,profile.(hvac.fans.schedule{f}));
                case 'Humidifier:Steam:Electric'
                    
                otherwise
                    disp('need to add components in manage_air_loops');
            end     
        end
        supply_air.T(i) = intermediate_air.T; supply_air.w(i) = intermediate_air.w; supply_air.m_dot(i) = intermediate_air.m_dot; supply_air.h(i) = intermediate_air.h;%update intermediate air as it passes between equipment
    end
end
end%Ends function sim_supply_equipment

function [zone_air,e_use] = sim_zone_equipment(building,supply_air,zone_flow,reheat_load,non_loop_mixed_air,T_air,T_supply_zone,e_use)
hvac = building.HVAC;
loops = length(hvac.zone_2_loop(:,1));
n = length(zone_flow);
zone_air.T = zeros(n,1);
zone_air.w = zeros(n,1);
zone_air.m_dot = zeros(n,1);
zone_air.h = zeros(n,1);
zone_air.m_dot = zone_flow;

for i = 1:1:loops%distribute through splitters to terminals (avoid divide by zero
    A = (hvac.zone_2_loop(i,:)>0);
    zone_air.T(A) = supply_air.T(i);
    zone_air.w(A) = supply_air.w(i);
    zone_air.h(A) = supply_air.h(i);
end
if ~isempty(building.HVAC.unitary_sys)
    for i = 1:1:length(building.HVAC.unitary_sys.name)
        if building.HVAC.unitary_sys.loop(i) == 0
            z = building.HVAC.unitary_sys.zone(i);
            zone_air.T(z) = non_loop_mixed_air.T(i);
            zone_air.w(z) = non_loop_mixed_air.w(i);
            zone_air.h(z) = non_loop_mixed_air.h(i);
            cp_supply = 1006 + 1860*zone_air.w(z);
            Q_treat = zone_air.m_dot(z).*cp_supply.*(T_supply_zone(z)-non_loop_mixed_air.T(i));
            load_heat = Q_treat.*(Q_treat>0);
            load_cool = Q_treat.*(Q_treat<0);
            air.T = zone_air.T(z);air.w = zone_air.w(z);air.m_dot = zone_air.m_dot(z);air.h = zone_air.h(z);%update intermediate air as it passes between equipment
            [air,e_use] = unitary_sys_operate(air,load_heat,load_cool,T_air,e_use,building.HVAC,i);
            zone_air.T(z) = air.T; zone_air.w(z) = air.w; zone_air.m_dot(z) = air.m_dot; zone_air.h(z) = air.h;
        end
    end
end

reheat_term = nonzeros((1:length(hvac.terminals.type))'.*strcmp(hvac.terminals.type,'Reheat'));
for i = 1:1:length(reheat_term)
    if reheat_load(i)>10%at least 10 W of re-heating
        z = hvac.terminals.zone(reheat_term(i));
        air.T = zone_air.T(z);air.w = zone_air.w(z);air.m_dot = zone_air.m_dot(z);air.h = zone_air.h(z);%update intermediate air as it passes between equipment
        name = hvac.terminals.reheat_coil_name{reheat_term(i)};
        switch hvac.terminals.reheat_coil_type{reheat_term(i)}
            case 'Coil:Heating:Fuel'
                [air,Q_coil] = heating_coil(hvac.coils.heating,name,'operation',air,reheat_load(i));
                e_use.heat_fuel = e_use.heat_fuel + Q_coil;
            case 'Coil:Heating:Electric'
                [air,Q_coil] = heating_coil(hvac.coils.heating,name,'operation',air,reheat_load(i));
                e_use.heat_elec = e_use.heat_elec + Q_coil;
            case 'Coil:Heating:Water'
                pl = building.plant.components.loop(nonzeros((1:length(building.plant.components.name))'.*strcmpi(name,building.plant.components.name)));%plant loop
                [air,water] = heating_coil(hvac.coils.heating,name,'operation',air,reheat_load(i),building.plant.loop.hot_water_supply_temperature(pl));
                e_use.heat_water.plant_loop(end+1) = pl;
                e_use.heat_water.return_flow(end+1) = water.m_dot;
                e_use.heat_water.return_temp(end+1) = water.T;
        end
        zone_air.T(z) = air.T; zone_air.w(z) = air.w; zone_air.m_dot(z) = air.m_dot; zone_air.h(z) = air.h;
    end
end
end%Ends function sim_zone_equipment

function [air,e_use] = unitary_sys_operate(air,load_heat,load_cool,T_air,e_use,hvac,i)
unit = hvac.unitary_sys;
if strcmp(unit.type{i},'heat_cool')
    name = unit.cool_coil{i};
    switch unit.cool_coil_type{i}
        case {'Coil:Cooling:DX:TwoSpeed';'Coil:Cooling:DX:SingleSpeed'}
            [air,P_coil] = cooling_coil(hvac.coils.cooling,name,'operation',hvac.curves,air,T_air,load_cool,air.w);
            e_use.cool_elec = e_use.cool_elec + P_coil;
        case {'Coil:Cooling:Water:DetailedGeometry';'Coil:Cooling:Water'}
            pl = building.plant.components.loop(nonzeros((1:length(building.plant.components.name))'.*strcmpi(name,building.plant.components.name)));%plant loop
            [air,water] = heating_coil(hvac.coils.heating,name,'operation',air,load_cool,building.plant.loop.cool_water_supply_temperature(pl));
            e_use.cool_water.plant_loop(end+1) = pl;
            e_use.cool_water.return_flow(end+1) = water.m_dot;
            e_use.cool_water.return_temp(end+1) = water.T;
    end
end
name = unit.heat_coil{i};
switch unit.heat_coil_type{i}
    case 'Coil:Heating:Fuel'
        [air,Q_coil] = heating_coil(hvac.coils.heating,name,'operation',air,load_heat);
        e_use.heat_fuel = e_use.heat_fuel + Q_coil;
    case 'Coil:Heating:Electric'
        [air,Q_coil] = heating_coil(hvac.coils.heating,name,'operation',air,load_heat);
        e_use.heat_elec = e_use.heat_elec + Q_coil;        
    case 'Coil:Heating:Water'
        pl = building.plant.components.loop(nonzeros((1:length(building.plant.components.name))'.*strcmpi(name,building.plant.components.name)));%plant loop
        [air,water] = heating_coil(hvac.coils.heating,name,'operation',air,load_heat,building.plant.loop.hot_water_supply_temperature(pl));
        e_use.heat_water.plant_loop(end+1) = pl;
        e_use.heat_water.return_flow(end+1) = water.m_dot;
        e_use.heat_water.return_temp(end+1) = water.T;
end
if strcmp(unit.type{i},'heat_cool') && ~isempty(unit.reheat_name{i})
    disp('need to add dehumid and re-heat capability to unitary systems')
end
end%Ends function unitary_sys_operate