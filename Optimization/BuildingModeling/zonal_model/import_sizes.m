function building = import_sizes(filename,building)
component = parse_eio(filename);
%% some component types could be combined, but I have kept everything seperate
for i = 1:1:length(component.name)
    switch component.type{i}
        case 'AirTerminal:SingleDuct:Uncontrolled'
            k = f_index(component.name{i},building.hvac.terminals.name);
            switch component.field{i}
                case 'Design Size Maximum Air Flow Rate [m3/s]'
                    building.hvac.terminals.max_flow(k) = component.value(i);
                case 'User-Specified Maximum Air Flow Rate [m3/s]'
                    if isnan(building.hvac.terminals.max_flow(k))%keep design value?
                        building.hvac.terminals.max_flow(k) =  component.value(i);
                    end
                otherwise
                    disp(strcat('need to load more parameters for component of type__',component.type{i}))
            end
        case 'AirTerminal:SingleDuct:VAV:Reheat'
            k = f_index(component.name{i},building.hvac.terminals.name);
            c = f_index(building.hvac.terminals.reheat_coil_name{k},building.coils.heating.name);
            switch component.field{i}
                case 'Design Size Maximum Air Flow Rate [m3/s]'
                    building.hvac.terminals.max_flow(k) = component.value(i);
                case 'Design Size Constant Minimum Air Flow Fraction'
                    building.hvac.terminals.min_frac(k) = component.value(i);
                case 'Design Size Minimum Air Flow Rate [m3/s]'
                    building.hvac.terminals.min_flow(k) = component.value(i);
                case 'Design Size Maximum Flow per Zone Floor Area during Reheat [m3/s-m2]'
                    building.hvac.terminals.max_reheat_flow_per_area(k) = component.value(i);
                case 'Design Size Maximum Flow Fraction during Reheat []'
                    building.hvac.terminals.max_reheat_frac(k) = component.value(i);
                case 'Design Size Maximum Reheat Water Flow Rate [m3/s]'
                    building.hvac.terminals.max_flow_water(k) = component.value(i);
                case 'Design Size Reheat Coil Sizing Air Volume Flow Rate [m3/s]'
                    if ~isfield(building.coils.heating,'rated_air_flow')
                        building.coils.heating.rated_air_flow = nan(length(building.coils.heating.name),1);
                    end
                    building.coils.heating.rated_air_flow(c,1) = component.value(i);
                case 'Design Size Reheat Coil Sizing Inlet Air Temperature [C]'
                    building.coils.heating.air_inlet_temperature(c) = component.value(i);
                case 'Design Size Reheat Coil Sizing Inlet Air Humidity Ratio [kgWater/kgDryAir]'
                    building.coils.heating.air_inlet_humidity_ratio(c,1) = component.value(i);
                case 'User-Specified Constant Minimum Air Flow Fraction'
                    building.hvac.terminals.min_frac(k) = component.value(i);
                otherwise
                    disp(strcat('need to load more parameters for component of type__',component.type{i}))
            end
        case 'Coil:Heating:Water'
            c = f_index(component.name{i},building.coils.heating.name);
            switch component.field{i}
                case 'Design Size Rated Capacity [W]'
                    building.coils.heating.design_capacity(c,1) = component.value(i);
                    cp_air = 1025; % J/kg*K
                    air_density = 1.204; %kg/m^3 %standard air density (1.204 kg/m3) adjusted for the local barometric pressure (standard barometric pressure corrected for altitude, ASHRAE 1997 HOF pg. 6.1).
                    building.coils.heating.rated_air_flow(c,1) = building.coils.heating.design_capacity(c,1)/(cp_air*air_density*(building.coils.heating.air_outlet_temperature(c,1)-building.coils.heating.air_inlet_temperature(c,1)));
                case 'Design Size Maximum Water Flow Rate [m3/s]'
                    building.coils.heating.rated_water_flow(c) = component.value(i);
                case 'Design Size U-Factor Times Area Value [W/K]'
                    building.coils.heating.UA(c) = component.value(i);
                case 'Nominal Total Capacity {W}'
                    building.coils.heating.capacity(c) = component.value(i);
                otherwise
                    disp(strcat('need to load more parameters for component of type__',component.type{i}))
            end
            if isnan(building.coils.heating.air_outlet_temperature(c))
                building.coils.heating.air_outlet_temperature(c) = building.hvac.design.Tsupply_h(al);
            end
        case {'Coil:Heating:Fuel';'Coil:Heating:Electric'}
            c = f_index(component.name{i},building.coils.heating.name);
            switch component.field{i}
                case 'Design Size Nominal Capacity [W]'
                    building.coils.heating.capacity(c,1) = component.value(i);
                otherwise
                    disp(strcat('need to load more parameters for component of type__',component.type{i}))
            end
        case 'Coil:Cooling:Water'
            c = f_index(component.name{i},building.coils.cooling.name);
            switch component.field{i}
                case 'Design Size Design Coil Load [W]'
                    building.coils.cooling.design_capacity(c,1) = component.value(i);
                case 'Design Size Design Water Flow Rate [m3/s]'
                    building.coils.cooling.rated_water_flow(c) = component.value(i);
                case 'Design Size Design Air Flow Rate [m3/s]'
                    building.coils.cooling.rated_air_flow(c) = component.value(i);
                case 'User-Specified Design Air Flow Rate [m3/s]'
                    if isnan(building.coils.cooling.rated_air_flow(c))%keep design value? override user input?
                        building.coils.cooling.rated_air_flow(c) = component.value(i);
                    end
                case 'Design Size Design Inlet Air Temperature [C]'
                    building.coils.cooling.air_inlet_temperature(c) = component.value(i);
                case 'Design Size Design Inlet Water Temperature [C]'
                    building.coils.cooling.water_inlet_temperature(c) = component.value(i);
                case 'Design Size Design Outlet Air Temperature [C]'
                    building.coils.cooling.air_outlet_temperature(c) = component.value(i);
                case 'Design Size Design Inlet Air Humidity Ratio'
                    building.coils.cooling.air_inlet_humididty(c) = component.value(i);
                case 'Design Size Design Outlet Air Humidity Ratio'
                    building.coils.cooling.air_outlet_humididty(c) = component.value(i);
                case 'Nominal Total Capacity {W}'
                    building.coils.cooling.capacity(c,1) = component.value(i);
                case 'Nominal Sensible Capacity {W}'
                    building.coils.cooling.sensible_capacity(c,1) = component.value(i);
                case 'Nominal Latent Capacity {W}'
                    building.coils.cooling.latent_capacity(c,1) = component.value(i);
                case 'Nominal Sensible Heat Ratio'
                    building.coils.cooling.sensible_heat_ratio(c,1) = component.value(i);
                case 'Nominal Coil UA Value {W/C}'
                    building.coils.cooling.UA(c,1) = component.value(i);
                case 'Nominal Coil Surface Area {m2}'
                    building.coils.cooling.area(c,1) = component.value(i);
                otherwise
                    disp(strcat('need to load more parameters for component of type__',component.type{i}))
            end
        case 'Coil:Cooling:Water:DetailedGeometry'
            c = f_index(component.name{i},building.coils.cooling.name);
            switch component.field{i}
                case 'Design Size Design Coil Load [W]'
                    building.coils.cooling.design_capacity(c,1) = component.value(i);
                case 'Design Size Maximum Water Flow Rate [m3/s]'
                    building.coils.cooling.max_flow(c) = component.value(i);
                case 'Design Size Number of Tubes per Row'
                    building.coils.cooling.tubes_per_row(c) = component.value(i);
                case 'Design Size Fin Diameter [m]'
                    building.coils.cooling.fin_diameter(c) = component.value(i);
                case 'Design Size Minimum Airflow Area [m2]'
                    building.coils.cooling.min_airflow_area(c) = component.value(i);
                case 'Design Size Fin Surface Area [m2]'
                    building.coils.cooling.fin_area(c) = component.value(i);
                case 'Design Size Total Tube Inside Area [m2]'
                    building.coils.cooling.area_inside(c) = component.value(i);
                case 'Design Size Tube Outside Surface Area [m2]'
                    building.coils.cooling.area_outside(c) = component.value(i);
                case 'Design Size Coil Depth [m]'
                    building.coils.cooling.coil_depth(c) = component.value(i);
                case 'Nominal Total Capacity {W}'
                    building.coils.cooling.capacity(c,1) = component.value(i);
                case 'Nominal Sensible Capacity {W}'
                    building.coils.cooling.sensible_capacity(c,1) = component.value(i);
                case 'Nominal Latent Capacity {W}'
                    building.coils.cooling.latent_capacity(c,1) = component.value(i);
                case 'Nominal Sensible Heat Ratio'
                    building.coils.cooling.sensible_heat_ratio(c,1) = component.value(i);
                otherwise
                    disp(strcat('need to load more parameters for component of type__',component.type{i}))
            end
        case 'Coil:Cooling:DX:TwoSpeed'
            c = f_index(component.name{i},building.coils.cooling.name);
            switch component.field{i}
                case 'Design Size High Speed Rated Air Flow Rate [m3/s]'
                    building.coils.cooling.rated_air_flow(c,1) = component.value(i);
                case 'Design Size High Speed Gross Rated Total Cooling Capacity [W]'
                    building.coils.cooling.rated_capacity(c) = component.value(i);
                case 'Design Size High Speed Rated Sensible Heat Ratio'
                    building.coils.cooling.rated_sensible_heat_ratio(c) = component.value(i);
                case 'Design Size Low Speed Rated Air Flow Rate [m3/s]'
                    building.coils.cooling.rated_air_flow2(c) = component.value(i);
                case 'Design Size Low Speed Gross Rated Total Cooling Capacity [W]'
                    building.coils.cooling.rated_capacity2(c) = component.value(i);
                case 'Design Size Low Speed Gross Rated Sensible Heat Ratio'
                    building.coils.cooling.rated_sensible_heat_ratio2(c) = component.value(i);
                otherwise
                    disp(strcat('need to load more parameters for component of type__',component.type{i}))
            end
        case 'Coil:Cooling:DX:SingleSpeed'
            c = f_index(component.name{i},building.coils.cooling.name);
            switch component.field{i}
                case 'Design Size Rated Air Flow Rate [m3/s]'
                    building.coils.cooling.rated_air_flow(c,1) = component.value(i);
                case {'Design Size Rated Total Cooling Capacity [W]';'Design Size Gross Rated Total Cooling Capacity [W]'}
                    building.coils.cooling.rated_capacity(c) = component.value(i);
                case {'Design Size Rated Sensible Heat Ratio';'Design Size Gross Rated Sensible Heat Ratio'}
                    building.coils.cooling.rated_sensible_heat_ratio(c) = component.value(i);
                otherwise
                    disp(strcat('need to load more parameters for component of type__',component.type{i}))
            end
        case 'AirLoopHVAC'
            al = f_index(component.name{i},building.hvac.loop.name);% air loop
            switch component.field{i}
                case 'Design Supply Air Flow Rate [m3/s]'
                    building.hvac.design.max_loop_flow(al,1) = component.value(i);
                otherwise
                    disp(strcat('need to load more parameters for component of type__',component.type{i}))
            end
        case 'AirLoopHVAC:UnitaryHeatCool'
            u = f_index(component.name{i},building.unitary_heat_cool.name);
            switch component.field{i}
                case 'Supply Air Flow Rate [m3/s]'
                    building.unitary_heat_cool.max_air_flow(u,1) = component.value(i);
                case 'Supply Air Flow Rate During Heating Operation [m3/s]'
                    building.unitary_heat_cool.heating_air_flow(u) = component.value(i);
                case 'Supply Air Flow Rate During Cooling Operation [m3/s]'
                    building.unitary_heat_cool.cooling_air_flow(u) = component.value(i);
                case 'Supply Air Flow Rate When No Cooling or Heating is Needed [m3/s]'
                    building.unitary_heat_cool.no_load_air_flow(u) = component.value(i);
                case 'Nominal Heating Capacity [W]'
                    building.unitary_heat_cool.capacity_heating(u,1) = component.value(i);
                case 'Nominal Cooling Capacity [W]'
                    building.unitary_heat_cool.capacity_cooling(u,1) = component.value(i);
                case 'Maximum Supply Air Temperature from Supplemental Heater [C]'
                    building.unitary_heat_cool.max_temperature(u) = component.value(i);
                case 'Fraction of Supply Air Flow That Goes Through the Controlling Zone'
                    building.unitary_heat_cool.frac_2_zone(u,1) = component.value(i);
                otherwise
                    disp(strcat('need to load more parameters for component of type__',component.type{i}))
            end
        case 'Controller:OutdoorAir'
            oa = f_index(component.name{i},building.hvac.outdoor_air.control);
            switch component.field{i}
                case 'Maximum Outdoor Air Flow Rate [m3/s]'
                    building.hvac.outdoor_air.max_flow(oa) = component.value(i);
                case 'Minimum Outdoor Air Flow Rate [m3/s]'
                    building.hvac.outdoor_air.min_flow(oa) = component.value(i);
                otherwise
                    disp(strcat('need to load more parameters for component of type__',component.type{i}))
            end
        case {'Fan:ConstantVolume';'Fan:VariableVolume';'Fan:OnOff'}
            f = f_index(component.name{i},building.fans.name);
            switch component.field{i}
                case 'Design Size Maximum Flow Rate [m3/s]'
                    building.fans.flow_rate(f) = component.value(i);
                otherwise
                    disp(strcat('need to load more parameters for component of type__',component.type{i}))
            end
        case {'PlantLoop';'CondenserLoop'}
            pl = f_index(component.name{i},building.plant_loop.name);%plant loop
            switch component.field{i}
                case 'Maximum Loop Flow Rate [m3/s]'
                    building.plant_loop.max_flow_rate(pl,1) = component.value(i);
                case {'Plant Loop Volume [m3]';'Condenser Loop Volume [m3]'}
                    building.plant_loop.loop_volume(pl,1) = component.value(i);
                otherwise
                    disp(strcat('need to load more parameters for component of type__',component.type{i}))
            end
        case {'Pump:VariableSpeed';'Pump:ConstantSpeed'}
            p = f_index(component.name{i},building.pump.name);
            switch component.field{i}
                case 'Design Flow Rate [m3/s]'
                    building.pump.design_flow(p) = component.value(i);
                case 'Design Power Consumption [W]'
                    building.pump.design_power(p) = component.value(i);
                otherwise
                    disp(strcat('need to load more parameters for component of type__',component.type{i}))
            end
        case 'WaterHeater:Mixed'
            wh = f_index(component.name{i},building.water_heater.name);
            switch component.field{i}
                case 'Use Side Design Flow Rate [m3/s]'
                    building.water_heater.use_side_design_flow(wh) = component.value(i);
                case 'Source Side Design Flow Rate [m3/s]'
                    building.water_heater.source_side_design_flow(wh) = component.value(i);
                otherwise
                    disp(strcat('need to load more parameters for component of type__',component.type{i}))
            end
        case 'Controller:WaterCoil'
            c = f_index(component.name{i},building.controller.name);
            switch component.field{i}
                case 'Maximum Actuated Flow [m3/s]'
                    building.controller.max_actuated_flow(c) = component.value(i);
                case 'Controller Convergence Tolerance'
                    building.controller.tolerance(c) = component.value(i);
                otherwise
                    disp(strcat('need to load more parameters for component of type__',component.type{i}))
            end
        case {'Chiller:Electric:ReformulatedEIR';'Chiller:Electric:EIR'}
            c = f_index(component.name{i},building.chiller.name);
            switch component.field{i}
                case 'Design Size Reference Chilled Water Flow Rate [m3/s]'
                    building.chiller.chilled_water_flow_ref(c) = component.value(i);
                case 'Design Size Reference Capacity [W]'
                    building.chiller.capacity(c) = component.value(i);
                case 'Design Size Reference Condenser Water Flow Rate [m3/s]'
                    building.chiller.condenser_water_flow_ref(c) = component.value(i);
                case 'Design Size Reference Condenser Fluid Flow Rate  [m3/s]'
                    building.chiller.condenser_water_flow_ref(c) = component.value(i);
                otherwise
                    disp(strcat('need to load more parameters for component of type__',component.type{i}))
            end         
        case 'Boiler:HotWater'
            b = f_index(component.name{i},building.boiler.name);
            switch component.field{i}
                case 'Design Size Nominal Capacity [W]'
                    building.boiler.capacity(b) = component.value(i);
                case 'Design Size Design Water Flow Rate [m3/s]'
                    building.boiler.design_flow_rate(b) = component.value(i);
                otherwise
                    disp(strcat('need to load more parameters for component of type__',component.type{i}))
            end
        case {'CoolingTower:SingleSpeed';'CoolingTower:TwoSpeed';'CoolingTower:VariableSpeed:Merkel'}
            t = f_index(component.name{i},building.cooling_tower.name);
            switch component.field{i}
                case 'Design Water Flow Rate [m3/s]'
                    building.cooling_tower.water_flow(t) = component.value(i);
                case 'Fan Power at Design Air Flow Rate [W]'
                    building.cooling_tower.fan_power(t) = component.value(i);
                case 'Design Air Flow Rate [m3/s]'
                    building.cooling_tower.air_flow(t) = component.value(i);
                case 'U-Factor Times Area Value at Design Air Flow Rate [W/C]'
                    building.cooling_tower.UA(t) = component.value(i);
                case 'Free Convection Regime Air Flow Rate [m3/s]'
                    building.cooling_tower.free_convect_flow(t) = component.value(i);
                case 'Free Convection U-Factor Times Area Value [W/K]'
                    building.cooling_tower.free_convect_UA(t) = component.value(i);
                otherwise
                    disp(strcat('need to load more parameters for component of type__',component.type{i}))
            end
        case 'Humidifier:Steam:Electric'
            %%need to decide to keep design or user specified sizes
            h = f_index(component.name{i},building.humidifiers.name);
            switch component.field{i}
                case 'Design Size Nominal Capacity Volume [m3/s]'
                    building.humidifiers.flow_capacity(h) = component.value(i);
                case 'User-Specified Nominal Capacity Volume [m3/s]'
%                     building.humidifiers.flow_capacity(h) = component.value(i);
                case 'Design Size Rated Power [W]'
                    building.humidifiers.max_power(h) = component.value(i);
                case 'User-Specified Rated Power [W]'
%                     building.humidifiers.max_power(h) = component.value(i);
                otherwise
                    disp(strcat('need to load more parameters for component of type__',component.type{i}))
            end
        case 'ZoneHVAC:UnitHeater'
            h = f_index(component.name{i},building.unitary_sys.name);
            switch component.field{i}
                case 'Design Size Maximum Supply Air Flow Rate [m3/s]'
                    building.unitary_sys.max_air_flow(h) = component.value(i);
                otherwise
                    disp(strcat('need to load more parameters for component of type__',component.type{i}))
            end
        case 'ZoneHVAC:FourPipeFanCoil'
            f = f_index(component.name{i},building.unitary_sys.name);
            switch component.field{i}
                case 'Design Size Maximum Supply Air Flow Rate [m3/s]'
                    building.unitary_sys.max_air_flow(f) = component.value(i);
                case 'Design Size Maximum Hot Water Flow [m3/s]'
                    building.unitary_sys.max_cold_water(f) = component.value(i);
                case 'Design Size Maximum Cold Water Flow [m3/s]'
                    building.unitary_sys.max_hot_water(f) = component.value(i);
                otherwise
                    disp(strcat('need to load more parameters for component of type__',component.type{i}))
            end
        case 'ZoneHVAC:PackagedTerminalAirConditioner'
            p = f_index(component.name{i},building.unitary_sys.name);
            switch component.field{i}
                case 'Design Size Cooling Supply Air Flow Rate [m3/s]'
                    building.unitary_sys.cooling_air_flow(p) = component.value(i);
                case 'Design Size Heating Supply Air Flow Rate [m3/s]'
                    building.unitary_sys.heating_air_flow(p) = component.value(i);
                case 'Design Size No Load Supply Air Flow Rate [m3/s]'
                    building.unitary_sys.no_load_air_flow(p) = component.value(i);
                case 'Supply Air Flow Rate When No Cooling or Heating is Needed [m3/s]'
                    building.unitary_sys.no_load_air_flow(p) = component.value(i);
                case 'Design Size Outdoor Air Flow Rate During Cooling Operation [m3/s]'
                    building.unitary_sys.cooling_outdoor_air_flow(p) = component.value(i);
                case 'Design Size Outdoor Air Flow Rate During Heating Operation [m3/s]'
                    building.unitary_sys.heating_outdoor_air_flow(p) = component.value(i);
                case 'Design Size Outdoor Air Flow Rate When No Cooling or Heating is Needed [m3/s]'
                    building.unitary_sys.no_load_outdoor_air_flow(p) = component.value(i);
                otherwise
                    disp(strcat('need to load more parameters for component of type__',component.type{i}))
            end
        otherwise
            disp(strcat('need to load sizes for component of type__',component.type{i}))
    end
end

[building.coils.cooling.water_inlet_temperature,building.coils.cooling.water_outlet_temperature] = water_setpoint(building,'cooling');
[building.coils.heating.water_inlet_temperature,building.coils.heating.water_outlet_temperature] = water_setpoint(building,'heating');
end%Ends import_sizes

function component = parse_eio(filename)
[cat,heading,z_info] = read_eio(filename);
comp = f_index('Component Sizing Information',cat);
c_objects = z_info{comp};

for i = 1:1:length(c_objects(:,1))
    component.type(i,1) = c_objects(i,1);
    component.name(i,1) = c_objects(i,2);
    component.field(i,1) = c_objects(i,3);
    component.value(i,1) = str2double(c_objects{i,4});
end

coil_heat = f_index('Water Heating Coil Capacity Information',cat);
if ~isempty(coil_heat)
    hc_objects = z_info{coil_heat};
    for j = i+1:1:i+length(hc_objects(:,1))
        component.type(j,1) = hc_objects(j-i,1);
        component.name(j,1) = hc_objects(j-i,2);
        component.field(j,1) = heading{coil_heat}(3);
        component.value(j,1) = str2double(hc_objects{j-i,3});
    end
end

coil_cool = f_index('Water Cooling Coil Capacity Information',cat);
if ~isempty(coil_cool)
    cc_objects = z_info{coil_cool};
    f = length(heading{coil_cool})-2;
    k = j;
    for p = j+1:1:j+length(cc_objects(:,1))
        component.type(k+1:k+f,1) = cc_objects(p-j,1);
        component.name(k+1:k+f,1) = cc_objects(p-j,2);
        for n = 1:1:f
            component.field(k+n,1) = heading{coil_cool}(n+2);
            component.value(k+n,1) = str2double(cc_objects{p-j,n+2});
        end
        k = k+f;
    end
end
end%ends function parse_eio

function [cat,heading,z_info] = read_eio(filename)
if ~contains(filename,'.eio')
    filename = strcat(filename,'.eio');
end
text = fileread(filename);
n_line = strfind(text,char(10));
n_l = length(n_line);
lines = cell(n_l,1);
s = 1;
for i = 1:1:length(n_line)
    lines(i) = {text(s:n_line(i)-2)};
    s = n_line(i)+1;
end

j = 0;
for i = 2:1:length(lines)-1
    com = strfind(lines{i},',');
    if strcmp(lines{i}(1),'!')
        j = j+1;
        cat(j,1) = {lines{i}(4:com(1)-2)};
        if length(com)==1
            headings = {rmv_spaces(lines{i}(com(1)+1:end))};
        else
            headings = {rmv_spaces(lines{i}(com(1)+1:com(2)-1))};
            for ii = 2:1:length(com)-1
                headings = [headings;rmv_spaces(lines{i}(com(ii)+1:com(ii+1)-1))];
            end
            headings = [headings;rmv_spaces(lines{i}(com(end)+1:end))];
        end
        heading(j,1) = {headings};
        z_info(j,1) = {cell(0,length(com))};
    else
        cat_i = f_index(rmv_spaces(lines{i}(1:com(1)-1)),cat);
        if isempty(cat_i)
            cat_i = length(cat);
        end
        k = length(z_info{cat_i}(:,1))+1;
        if length(com)==1
            z_info{cat_i,1}{k,1} = rmv_spaces(lines{i}(com(1)+1:end));
        else
            for ii = 1:1:length(com)-1
                z_info{cat_i,1}{k,ii} = rmv_spaces(lines{i}(com(ii)+1:com(ii+1)-1));
            end
            z_info{cat_i,1}{k,length(com)} = rmv_spaces(lines{i}(com(end)+1:end));
        end
    end
end
end%Ends function read_eio

function [Tw_in,Tw_out] = water_setpoint(building,type)
Tw_out = [];
Tw_in = [];
if any(ismember(building.coils.(type).type,{'Water';'Water:DetailedGeometry'}))
    Tw_out = nan(length(building.coils.(type).name),1);
    Tw_in = Tw_out;
    if isfield(building.coils.(type),'water_inlet_temperature')
        Tw_in = building.coils.(type).water_inlet_temperature;
    end
    for c = 1:1:length(building.coils.(type).name)
        pl = 0;
        name = building.coils.(type).name{c};
        i = f_index(name,building.plant_demand_equip.name);
        if ~isempty(i)
            pl = building.plant_demand_equip.loop(i);
        else
            i = f_index(name,building.plant_supply_equip.name);
            if ~isempty(i)
                pl = building.plant_supply_equip.loop(i);
            end
        end
        switch building.coils.(type).type{c}
            case {'Water';'Water:DetailedGeometry'}
                if isnan(Tw_in(c))
                    Tw_in(c) = building.plant_loop.exit_temperature(pl);
                end
                switch type
                    case 'cooling'
                        Tw_out(c) = Tw_in(c) + building.plant_loop.temperature_difference(pl);
                    case 'heating'
                        Tw_out(c) = Tw_in(c) - building.plant_loop.temperature_difference(pl);
                end
        end
    end
end
end%Ends function water_setpoint

function name = rmv_spaces(name)
while ~isempty(name) && strcmp(name(1),' ')
    name = name(2:end);
end
while ~isempty(name) && strcmp(name(end),' ')
    name = name(1:end-1);
end
end%Ends function rmv_spaces