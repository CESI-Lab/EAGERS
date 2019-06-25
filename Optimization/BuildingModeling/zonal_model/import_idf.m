function building = import_idf(filename)
[objects,ground_temps,convection] = read_objects(filename);
%name
i = f_index('Building',objects.type);
building.name = objects.name{i};
prop = {'north_axis';'terrain';'load_tol';'temp_tol';'solar_distribution';'max_warm_up';'min_warm_up';};
building.site = load_object(objects,'Building',{''},{prop},0);
%ground reflectance
A = f_index('Site:GroundReflectance',objects.type);
if ~isempty(A)
    disp('building has ground reflectance example input')
end
building.site.ground_reflect = 0;

%date
i = f_index('RunPeriod',objects.type);
building.sim_date = import_date(objects.value{i},1);

%location
building.location = {};
i = f_index('Site:Location',objects.type);
building.location.Latitude = str2double(objects.value{i}{1});
building.location.Longitude = str2double(objects.value{i}{2});
building.location.TimeZone = str2double(objects.value{i}{3});
building.location.Elevation = str2double(objects.value{i}{4});

%convection algorithm
building.convection = convection;

%schedules
A = f_index('Schedule:Compact',objects.type);
building.schedule = import_schedules(objects.name(A),objects.value(A));

%holidays
A = f_index('RunPeriodControl:SpecialDays',objects.type);
building.holidays.name = objects.name(A);
for k = 1:1:length(A)
    building.holidays.start(k,1) = objects.value{A(k)}(1);
    building.holidays.duration(k,1) = str2double(objects.value{A(k)}{2});
    building.holidays.type(k,1) = objects.value{A(k)}(3);
end


%Materials and construction
prop1 = {'roughness';'thickness';'conductivity';'density';'specheat';'thermal_absorptance';'solar_absorptance';'visible_absorptance';};
prop2 = {'roughness';'thermal_resistance';'thermal_absorptance';'solar_absorptance';'visible_absorptance';};
building.material = load_object(objects,'Material',{'';'NoMass';},{prop1,prop2},0);

prop1 = {'u_factor';'solar_heat_gain';'visible_transmittance';};
prop2 = {'optical_data_type';'spectral_data_set';'thickness';'solar_transmittance';'reflectance_front';'reflectance_back';'transmittance_visible';'reflectance_visible_front';'reflectance_visible_back';'transmittance_infrared';'emittance_front';'emittance_back';'thermal_conductivity';'dirt_correction';'solar_diffusing';};
prop3 = {'gas';'gap_thickness';};
building.window_material = load_object(objects,'WindowMaterial',{'SimpleGlazingSystem';'Glazing';'Gas'},{prop1,prop2,prop3},0);

A = f_index('Construction',objects.type);
i_nodes = 4;
[construction,window_construction] = read_construction(objects.name(A),objects.value(A),building.material,building.window_material,i_nodes);

%surfaces
building.surfaces = import_idf_surfaces(objects,construction);
building.windows = import_window(objects,window_construction);
building.doors = import_door(objects,construction);

%%frames
prop1 = {'width';'projection_out';'projection_in';'conductance';'conductance_ratio_edge_center';'solar_absoprtance';'visible_absorptance';'emissivity';'divider_type';'divider_width';'horizontal_dividers';'vertical_dividers';'divider_projection_out';'divider_projection_in';'divider_conductance';'divider_conductance_ratio';'divider_solar_absoprtnace';'divider_visible_absorptance';'divider_emissivity'};
building.window_frames = load_object(objects,'WindowProperty',{'FrameAndDivider'},{prop1},0);

%zone parameters & geometry of surfaces
prop1 = {'relative_north';'x_origin';'y_origin';'z_origin';'type';'multiplier';'ceil_height';'volume';'floor_area';'inside_convection';'outside_convection';'incl_in_total_area'};
building.zones = load_object(objects,'',{'Zone'},{prop1},0);

[building.zones,building.surfaces,building.windows] = zone_geometry(building.zones,building.surfaces,building.windows);

%infiltration
building.infiltration = import_infiltration(objects,building.zones);

%mixing
building.mixing = import_mixing(objects,building.zones,building.schedule);

%%occupancy
A = f_index('People',objects.type);
building.occupancy = import_occupancy(objects.name(A),objects.value(A),building.zones);
building.zones.max_occupancy = zeros(length(building.zones.name),1);
for i = 1:1:length(building.occupancy.zone)
    z_i = building.occupancy.zone(i);
    building.zones.max_occupancy(z_i) = building.zones.max_occupancy(z_i) + building.occupancy.nominal(i)*building.zones.floor_area(z_i); 
end

%%lighting and equipment in zones
A = f_index('Lights',objects.type);
building.lighting_internal = import_lighting(objects.name(A),objects.value(A),building.zones,building.occupancy);
A = f_index('ElectricEquipment',objects.type);
building.plug_load = import_equipment(objects.name(A),objects.value(A),building.zones,building.occupancy);
A = f_index('GasEquipment',objects.type);
building.gas_load = import_equipment(objects.name(A),objects.value(A),building.zones,building.occupancy);
%%external lighting and equipment
prop1 = {'schedule';'nominal';'control';'category';};
prop2 = {'fuel';'schedule';'nominal';'category';};
building.exterior = load_object(objects,'Exterior',{'Lights';'FuelEquipment';},{prop1,prop2},0);

%refrigeration cases & compressors 
prop = {'schedule';'zone';'rated_ambient_T';'rated_ambient_RH';'capacity_per_length';'latent_heat_ratio';'runtime_fraction';'length';'temperature';'latent_credit_curve_type';'latent_credit_curve_name';'fan_power_per_length';'operating_fan_power_per_length';'standard_lighting_per_unit_length';'installed_standard_lighting_per_unit_length';'light_schedule';'light_to_case';'anti_sweat_heater_per_length';'minimum_anti_sweat_per_length';'anti_sweat_control';'humidity_at_zero_percent';'height';'anti_sweat_heat_to_case';'defrost_power_per_length';'defrost_type';'defrost_schedule';'defrost_drip_down_schedule';'defrost_curve_type';'defrost_curve_name';'return_air_frac';'restock_schedule';'case_credit_fraction_schedule'};
building.cases = load_object(objects,'Refrigeration',{'Case'},{prop},0);
building.walk_ins = [];%need to add this object type
prop = {'heat_reject_location';'design_COP';'COP_curve';'fan_power';'fan_power_temperature_curve';'condensor_type';'water_inlet';'water_outlet';'water_loop_type';'water_condensor_temperature_schedule';'water_flow_rate';'water_max_flow';'water_max_outlet_T';'water_min_inlet_T';'evap_condenser_schedule';'evap_condenser_effectiveness';'evap_condenser_air_flow';'basin_heater_capacity';'basin_heater_setpoint';'water_pump_power';'water_supply_tank';'condenser_air_inlet';'end_use';'case_name';};
building.racks = load_object(objects,'Refrigeration',{'CompressorRack'},{prop},0);
D = f_index('Refrigeration:CaseAndWalkInList',objects.type);
building.racks = rack_cases(building.racks,building.cases,objects.name(D),objects.value(D));

%water use
prop1 = {'end_use_category';'peak_flow';'flow_schedule';'target_temperature_schedule';'hot_supply_temperature_schedule';'cold_supply_temperature_schedule';'zone';'sensible_frac_schedule';'latent_frac_schedule'};
prop2 = {'inlet';'outlet';'supply_tank';'reclaimation_tank';'hot_schedule';'cold_schedule';'drain_HX_type';'drain_HX_destination';'drain_HX_UA';'equip_name';};
building.water_use = load_object(objects,'WaterUse',{'Equipment';'Connections'},{prop1,prop2},0);

%HVAC systems
[building.hvac,e_list] = import_hvac(objects,building.zones.name,building.schedule);
[building.air_demand_nodes,building.air_demand_equip] = setup_nodes(building.hvac.components,building.hvac.loop,'d');
[building.air_supply_nodes,building.air_supply_equip] = setup_nodes(building.hvac.components,building.hvac.loop,'s');
for i = 1:1:length(building.hvac.loop.name)
    b = f_index(i,building.hvac.branches.loop);
    components = building.hvac.branches.component_name{b};
    for j = 1:1:length(components)
        k = f_index(components{j},building.air_supply_equip.name);
        building.air_supply_equip.branch_order(k,1) = j;
    end
end
    
%Plant loops
[building.plant_loop,components] = import_plant(objects);
[building.plant_demand_nodes,building.plant_demand_equip] = setup_nodes(components,building.plant_loop,'d');%% Demand Side
[building.plant_supply_nodes,building.plant_supply_equip] = setup_nodes(components,building.plant_loop,'s');%% Supply side

%pumps, chillers, boilers, humidifiers, cooling tower, water heater
prop1 = {'inlet';'outlet';'design_flow';'design_head';'design_power';'motor_eff';'fluid_ineff_frac';'c1';'c2';'c3';'c4';'design_min_flow';'control_type';'schedule';};
prop2 = {'inlet';'outlet';'design_flow';'design_head';'design_power';'motor_eff';'fluid_ineff_frac';'control_type';'schedule';};
building.pump = load_object(objects,'Pump',{'VariableSpeed';'ConstantSpeed'},{prop1,prop2},0);

prop1 = {'capacity';'COP_ref';'chilled_water_temperature_ref';'condenser_water_temperature_ref';'chilled_water_flow_ref';'condenser_water_flow_ref';'capacity_temperature_curve';'electric_in_out_temperature_curve';'electric_in_out_part_load_curve';'min_part_load';'max_part_load';'opt_part_load';'min_unload_ratio';'chilled_water_inlet';'chilled_water_outlet';'condenser_inlet';'condenser_outlet';'condenser_type';'fan_power_ratio';'frac_power_2_condenser';'min_chilled_water_temperature';'flow_mode';'heat_recovery_flow';'heat_recovery_inlet';'heat_recovery_outlet';'sizing_factor';};
prop2 = {'capacity';'COP_ref';'chilled_water_temperature_ref';'condenser_water_temperature_ref';'chilled_water_flow_ref';'condenser_water_flow_ref';'capacity_temperature_curve';'electric_in_out_temperature_curve';'electric_in_out_curve_type';'electric_in_out_part_load_curve';'min_part_load';'max_part_load';'opt_part_load';'min_unload_ratio';'chilled_water_inlet';'chilled_water_outlet';'condenser_inlet';'condenser_outlet';'frac_power_2_condenser';'min_chilled_water_temperature';'flow_mode';'heat_recovery_flow';'heat_recovery_inlet';'heat_recovery_outlet';'sizing_factor';};
building.chiller = load_object(objects,'Chiller',{'Electric:EIR';'Electric:ReformulatedEIR';'ConstantCOP';'Absorption';'Absorption:Indirect';'CombustionTurbine';'EngineDriven';},{prop1,prop2},0);%need examples of types 3-10 to create property lists

prop1 = {'fuel';'capacity';'efficiency';'eff_temperature_curve';'normalized_eff_curve';'design_water_temperature';'design_flow_rate';'min_part_load';'max_part_load';'optimal_part_load';'water_inlet';'water_outlet';'max_water_temperature';'flow_mode';'parasitic_electric_load';'sizing_factor';};
building.boiler = load_object(objects,'Boiler',{'HotWater'},{prop1},0);

prop1 = {'schedule';'flow_capacity';'max_power';'fan_power';'standby_power';};
building.humidifiers = load_object(objects,'Humidifier',{'Steam:Electric'},{prop1},0);

prop1 = {'schedule';'inlet';'outlet';'fan_schedule';'max_temperature';'cooling_air_flow';'heating_air_flow';'no_load_air_flow';'thermostat';'fan_type';'fan_name';'fan_location';'heat_coil_type';'heat_coil';'cool_coil_type';'cool_coil';'dehumid_control';'reheat_type';'reheat_name';};
building.unitary_heat_cool = load_object(objects,'AirLoopHVAC',{'UnitaryHeatCool';},{prop1},0);%'UnitaryHeatPump:AirToAir';
for i = 1:1:length(building.unitary_heat_cool.name)
    building.unitary_heat_cool.loop(i,1) = f_index(building.unitary_heat_cool.outlet{i},building.hvac.loop.supply_outlet);
end
        
prop1 = {'water_inlet';'water_outlet';'water_flow';'air_flow';'fan_power';'UA';'free_convect_flow';'free_convect_factor';'free_convect_UA';'free_convect_UA_factor';'performance_method';'heat_reject_capacity';'nominal_capacity';'free_convect_capacity';'free_convect_size_factor';'des_inlet_dry_bulb';'des_inlet_wet_bulb';'des_approach_temp';'des_range_temp';'basin_heater_capacity';'basin_heater_setpoint';'basin_heater_schedule';'evap_loss_mode';'evap_loss_factor';'drift_loss';'blowdown_mode';'blowdown_ratio';'blowdown_makeup_schedule';'supply_water_tank';'outdoor_air_inlet';'capacity_control';'number_of_cells';'cell_control';'cell_min_water_frac';'cell_max_water_frac';'sizing_factor';};
building.cooling_tower = load_object(objects,'CoolingTower',{'SingleSpeed';'TwoSpeed';'VariableSpeed:Merkel'},{prop1},0);%need examples of these other cooling tower types

prop1 = {'tank_volume';'temperature_schedule';'deadband_temperature_dif';'max_temperature';'heater_control_type';'heater_maximum_capacity';'heater_minimum_capacity';'heater_ignition_minimum_flow';'heater_ignition_delay';'heater_fuel';'heater_efficiency';'part_load_curve';'off_cycle_fuel_use';'off_cycle_fuel_type';'off_cycle_frac2tank';'on_cycle_fuel_use';'on_cycle_fuel_type';'on_cycle_frac2tank';'ambient_temp_indicator';'ambient_temp_schedule';'ambient_temp_zone';'ambient_temp_outdoor_node';'off_cycle_loss2ambient';'off_cycle_loss2zone';'on_cycle_loss2ambient';'on_cycle_loss2zone';'peak_use_flow';'use_flow_schedule';'cold_water_supply_temp_schedule';'use_side_inlet';'use_side_outlet';'use_side_effectiveness';'source_side_inlet';'source_side_outlet';'source_side_effectiveness';'use_side_design_flow';'source_side_design_flow';'indirect_water_heating_recovery_time';};
building.water_heater = load_object(objects,'WaterHeater',{'Mixed'},{prop1},0);

prop1 = {'control_variable';'action';'actuator_variable';'sensor_node';'actuator_node';'tolerance';'max_actuated_flow';'min_actuated_flow';};
building.controller = load_object(objects,'Controller',{'WaterCoil'},{prop1},0);

%zone/plant managers
prop1 = {'control_variable';'schedule';'node';};
prop2 = {'control_variable';'reference_setpoint_node';'fan_inlet_node';'fan_outlet_node';'node';};
prop3 = {'control_variable';'min_temp';'max_temp';'zone';'zone_node';'zone_inlet_node';'node'};
prop4 = {'node';'control_zone_node';};
prop6 = {'control_variable';'ref_temp_type';'offset_temperature';'max_temperature';'min_temperature';'node'};
prop7 = {'control_variable';'setpoint_at_low_temp';'low_temp';'setpoint_at_high_temp';'high_temp';'node'};
building.manager = load_object(objects,'SetpointManager',{'Scheduled';'MixedAir';'SingleZone:Reheat';'SingleZone:Humidity:Minimum';'SingleZone:Humidity:Maximum';'FollowOutdoorAirTemperature';'OutdoorAirReset';},{prop1,prop2,prop3,prop4,prop4,prop6,prop7},0);
%replace node_lists with nodes
n_list = load_object(objects,'NodeList',{''},{{'node'}},1);
for i = 1:1:length(building.manager.name)
    j = f_index(building.manager.node{i},n_list.name);
    if ~isempty(j)
        building.manager.node{i} = n_list.node{j};
    end
end
prop1 = {'zone';'control_type_schedule';'control_type';'control_name';};
prop2 = {'zone';'humidify_schedule';'dehumidify_schedule'};
building.zone_controls = load_object(objects,'ZoneControl',{'Thermostat';'Humidistat';},{prop1,prop2},2);
building.thermostat = load_object(objects,'ThermostatSetpoint',{'SingleHeating';'SingleCooling';'SingleHeatingOrCooling';'DualSetpoint'},{{'heating'},{'cooling'},{'heating'},{'heating';'cooling'}},0);

%zone supply flow temperature and humidity setpoints
A = f_index('DesignSpecification:OutdoorAir',objects.type);
outdoor_air = import_idf_outdoor_air(objects.name(A),objects.value(A));
prop1 = {'cooling_T_method';'cooling_T_des';'cooling_dT_des';'heating_T_method';'heating_T_des';'heating_dT_des';'cooling_w_des';'heating_w_des';'outdoor_air_object';'heating_sizing';'cooling_sizing';'cooling_flow_des_method';'cooling_flow_des';'cooling_min_flow_per_area';'cooling_min_flow';'cooling_min_frac';'heating_flow_des_method';'heating_flow_des';'heating_min_flow_per_area';'heating_min_flow';'heating_min_frac';};
building.setpoints = load_object(objects,'Sizing',{'Zone'},{prop1},0);
for i = 1:1:length(building.setpoints.name)
    z = f_index(building.setpoints.name{i},building.zones.name);
    building.setpoints.zone(i,1) = z;
    oa = f_index(building.setpoints.outdoor_air_object{i},outdoor_air.name);
    building.setpoints.outdoor_flow_method(i,1) = outdoor_air.method(oa);
    building.setpoints.outdoor_flow_value(i,1) = outdoor_air.value(oa);   
end

%Unitary Systems
prop1 = {'schedule';'inlet';'outlet';'fan_type';'fan_name';'max_air_flow';'heat_coil_type';'heat_coil';'fan_schedule';'fan_without_heating';'max_hot_water';'min_hot_water';'tolerance';'availability';};
prop2 = {'schedule';'inlet';'outlet';'outdoor_air_type';'outdoor_air';'cooling_air_flow';'heating_air_flow';'no_load_air_flow';'cooling_outdoor_air_flow';'heating_outdoor_air_flow';'no_load_outdoor_air_flow';'fan_type';'fan_name';'heat_coil_type';'heat_coil';'cool_coil_type';'cool_coil';'fan_placement';'fan_schedule';};
prop3 = {'schedule';'capacity_control';'max_air_flow';'low_flow_ratio';'medium_flow_ratio';'outdoor_air_flow';'outdoor_air_schedule';'inlet';'outlet';'outdoor_air_mixer_type';'outdoor_air_mixer';'fan_type';'fan_name';'cool_coil_type';'cool_coil';'max_cold_water';'min_cold_water';'cooling_tolerance';'heat_coil_type';'heat_coil';'max_hot_water';'min_hot_water';'heating_tolerance';};
building.unitary_sys = load_object(objects,'ZoneHVAC',{'UnitHeater';'PackagedTerminalAirConditioner';'FourPipeFanCoil';},{prop1,prop2,prop3},0);%'PackagedTerminalHeatPump';'WindowAirConditioner';
building.unitary_sys.zone = zeros(length(building.unitary_sys.name),1);
for i = 1:1:length(building.unitary_sys.name)
    k = f_index(building.unitary_sys.name{i},e_list.e_name);
    building.unitary_sys.zone(i,1) = e_list.zone(k);
end

%Fans
prop1 = {'schedule';'fan_efficiency';'pressure_rise';'flow_rate';'motor_efficiency';'motor_frac';'inlet';'outlet';'end_use';};
prop2 = {'schedule';'fan_efficiency';'pressure_rise';'flow_rate';'inlet';'outlet';'end_use';};
prop3 = {'schedule';'fan_efficiency';'pressure_rise';'flow_rate';'method';'min_flow_frac';'min_flow_rate';'motor_efficiency';'motor_frac';'a1';'a2';'a3';'a4';'a5';'inlet';'outlet';'end_use';};
prop4 = {'schedule';'fan_efficiency';'pressure_rise';'flow_rate';'motor_efficiency';'motor_frac';'inlet';'outlet';'power_speed_curve';'efficiency_speed_curve';'end_use';};
building.fans = load_object(objects,'Fan',{'ConstantVolume';'ZoneExhaust';'VariableVolume';'OnOff'},{prop1,prop2,prop3,prop4},0);
building.fans.exhaust = false(length(building.fans.name),1);
building.fans.exhaust_zone = zeros(length(building.fans.name),1);
for i = 1:1:length(building.fans.name)
    building.fans.schedule(i) = {fix_schedule_name(building.fans.schedule{i},building.schedule)};
    k = f_index(building.fans.name{i},e_list.e_name);
    if ~isempty(k)
        building.fans.exhaust(i) = true;
        building.fans.exhaust_zone(i) = e_list.zone(k);
    end
end
        
%Efficiency and other curves
prop1 = {'a0';'a1';'a2';'min_x';'max_x';};
prop2 = {'a0';'a1';'a2';'a3';'min_x';'max_x';};
prop3 = {'a0';'a1';'a2';'b1';'b2';'ab';'min_x';'max_x';'min_y';'max_y';};
prop4 = {'a0';'a1';'a2';'b1';'b2';'ab';'a3';'b3';'aab';'abb';'min_x';'max_x';'min_y';'max_y';};
prop5 = {'a0';'a1';'a2';'b1';'b2';'ab';'a3';'b3';'aab';'abb';'aabb';'cb3';'min_x';'max_x';'min_y';'max_y';};
building.curves = load_object(objects,'Curve',{'Quadratic';'Cubic';'Biquadratic';'Bicubic'},{prop1,prop2,prop3,prop4,prop5},0);

%cooling coils
prop1 = {'schedule';'rated_capacity';'rated_sensible_heat_ratio';'rated_COP';'rated_air_flow';'static_air_pressure';'inlet';'outlet';'capacity_v_temperature_curve';'capacity_v_flow_curve';'energy_input_v_temperature_curve';'energy_input_v_flow_curve';'part_load_curve';'rated_capacity2';'rated_sensible_heat_ratio2';'rated_COP2';'rated_air_flow2';'capacity_v_temperature_curve2';'energy_input_v_temperature_curve2';'condenser_inlet';};
prop2 = {'schedule';'rated_capacity';'rated_sensible_heat_ratio';'rated_COP';'rated_air_flow';'rated_fan_power_per_volume';'inlet';'outlet';'capacity_v_temperature_curve';'capacity_v_flow_curve';'energy_input_v_temperature_curve';'energy_input_v_flow_curve';'part_load_curve';'min_ambient_temp';};
prop3 = {'schedule';'rated_water_flow';'rated_air_flow';'water_inlet_temperature';'air_inlet_temperature';'air_outlet_temperature';'air_inlet_humididty';'air_outlet_humidity';'inlet_water';'outlet_water';'inlet';'outlet';'analysis_type';'heat_exchanger_configuration';};
prop4 = {'schedule';'max_flow';'area_outside';'area_inside';'fin_area';'min_airflow_area';'coil_depth';'fin_diameter';'fin_thickness';'tube_id';'tube_od';'tube_conductivity';'fin_conductivity';'fin_spacing';'tube_spacing';'tube_rows';'tubes_per_row';'inlet_water';'outlet_water';'inlet';'outlet';};
building.coils.cooling = load_object(objects,'Coil:Cooling',{'DX:TwoSpeed';'DX:SingleSpeed';'Water';'Water:DetailedGeometry'},{prop1,prop2,prop3,prop4},0);%'WaterToAirHeatPump:VariableSpeedEquationFit';

%heating coils
prop1 = {'schedule';'fuel';'efficiency';'capacity';'inlet';'outlet';'temperature_setpoint';};
prop2 = {'schedule';'UA';'rated_water_flow';'inlet_water';'outlet_water';'inlet';'outlet';'input_method';'capacity';'water_inlet_temperature';'air_inlet_temperature';'water_outlet_temperature';'air_outlet_temperature';'air_water_convection_ratio';};
prop3 = {'schedule';'efficiency';'capacity';'inlet';'outlet';};
building.coils.heating = load_object(objects,'Coil:Heating',{'Fuel';'Water';'Electric'},{prop1,prop2,prop3},0);

%ground temperatures
building.ground_temperatures = ground_temps;
%water mains
A = f_index('Site:WaterMainsTemperature',objects.type);
if ~isempty(A)
    building.water_main_temperature.method = objects.name(A);
    building.water_main_temperature.schedule = {fix_name(objects.value{A}{1})};
    building.water_main_temperature.outdoor_average = str2double(objects.value{A}{2});
    building.water_main_temperature.max_temp_difference = str2double(objects.value{A}{3});
else
    building.water_main_temperature = [];
end

A = f_index('EnvironmentalImpactFactors',objects.type);
if ~isempty(A)
    building.impact_factor.district_heat_efficiency = str2double(objects.name{A});
    building.impact_factor.district_cool_cop = str2double(objects.value{A}{1});
    building.impact_factor.steam_efficiency = str2double(objects.value{A}{2});
    building.impact_factor.carbon_factor_n2o = str2double(objects.value{A}{3});
    building.impact_factor.carbon_factor_ch4 = str2double(objects.value{A}{4});
    building.impact_factor.carbon_factor_co2 = str2double(objects.value{A}{5});
else
    building.impact_factor.steam_efficiency = 0.25;%default efficiency for "purchased heating" i.e. unconnected water equipment conversion to natural gas
end

%% sometimes the component list uses a coil system name instead of a coil name, make this edit now, so as to not worry about it later
prop1 = {'schedule';'inlet';'outlet';'sensor';'coil_type';'coil_name';};
coil_system = load_object(objects,'CoilSystem',{'Cooling:DX'},{prop1},0);
for i = 1:1:length(building.hvac.components.name)
    if any(strcmp(building.hvac.components.type{i},{'CoilSystem:Cooling:DX';'Coil:Cooling:Water:DetailedGeometry';'Coil:Cooling:Water'})) && ~any(strcmpi(building.hvac.components.name{i},building.coils.cooling.name))
        k = f_index(building.hvac.components.name{i},coil_system.name);
        building.hvac.components.name{i} = coil_system.coil_name{k};
    end
end
for i = 1:1:length(building.air_supply_equip.name)
    if any(strcmp(building.air_supply_equip.type{i},{'CoilSystem:Cooling:DX';'Coil:Cooling:Water:DetailedGeometry';'Coil:Cooling:Water'})) && ~any(strcmpi(building.air_supply_equip.name{i},building.coils.cooling.name))
        k = f_index(building.air_supply_equip.name{i},coil_system.name);
        building.air_supply_equip.name{i} = coil_system.coil_name{k};
    end
end

%% check if there are any components types missed
loaded = {'AirTerminal:SingleDuct:Uncontrolled';'AirTerminal:SingleDuct:VAV:Reheat';...
    'AirLoopHVAC';'AirLoopHVAC:ZoneSplitter';'AirLoopHVAC:ZoneMixer';'AirLoopHVAC:UnitaryHeatCool';...
    'AirLoopHVAC:ReturnPlenum';'AirLoopHVAC:SupplyPath';'AirLoopHVAC:ReturnPath';'Branch';'BranchList';...
    'AirLoopHVAC:OutdoorAirSystem';'AirLoopHVAC:ControllerList';...
    'AvailabilityManagerAssignmentList';'AvailabilityManager:NightCycle';...
    'AvailabilityManager:Scheduled';'Boiler:HotWater';'Building';'BuildingSurface:Detailed';...
    'Chiller:Electric:EIR';'Chiller:Electric:ReformulatedEIR';'Construction';...
    'Coil:Heating:Fuel';'Coil:Heating:Water';'Coil:Heating:Electric';...
    'Coil:Cooling:DX:TwoSpeed';'Coil:Cooling:DX:SingleSpeed';'Coil:Cooling:Water';...
    'Coil:Cooling:Water:DetailedGeometry';'CoilSystem:Cooling:DX';...
    'Connector:Splitter';'Connector:Mixer';'Sizing:Plant';'PlantLoop';'CondenserLoop';...
    'Controller:OutdoorAir';'Controller:WaterCoil';'CoolingTower:SingleSpeed';...
    'Curve:Quadratic';'Curve:Cubic';'Curve:Biquadratic';'Curve:Bicubic';...
    'DesignSpecification:OutdoorAir';'EnvironmentalImpactFactors';...
    'ElectricEquipment';'Exterior:Lights';'Exterior:FuelEquipment';...
    'Fan:ConstantVolume';'Fan:ZoneExhaust';'Fan:VariableVolume';'Fan:OnOff';...
    'FenestrationSurface:Detailed';'GasEquipment';'Humidifier:Steam:Electric';'InternalMass';...
    'Lights';'People';'Pump:VariableSpeed';'Pump:ConstantSpeed';...
    'Material';'Material:NoMass';'NodeList';...
    'Refrigeration:CaseAndWalkInList';'Refrigeration:CompressorRack';'Refrigeration:Case';...
    'RunPeriod';'RunPeriodControl:SpecialDays';...
    'Schedule:Compact';'SetpointManager:SingleZone:Reheat';'SetpointManager:FollowOutdoorAirTemperature';...
    'SetpointManager:MixedAir';'SetpointManager:SingleZone:Humidity:Maximum';...
    'SetpointManager:SingleZone:Humidity:Minimum';'SetpointManager:OutdoorAirReset';...
    'SetpointManager:Scheduled';'Site:Location';'Site:WaterMainsTemperature';...
    'Sizing:Zone';'Sizing:System';'SizingPeriod:DesignDay';'Sizing:Parameters';...
    'ThermostatSetpoint:DualSetpoint';...
    'WaterHeater:Mixed';'WaterUse:Equipment';'WaterUse:Connections';...
    'WindowMaterial:SimpleGlazingSystem';'WindowMaterial:Glazing';'WindowMaterial:Gas';'WindowProperty:FrameAndDivider';...
    'Zone';'ZoneHVAC:AirDistributionUnit';'ZoneHVAC:EquipmentList';'ZoneHVAC:EquipmentConnections';...
    'ZoneHVAC:UnitHeater';'ZoneHVAC:UnitaryHeatCool';'ZoneHVAC:PackagedTerminalAirConditioner';'ZoneHVAC:FourPipeFanCoil';...
    'ZoneInfiltration:DesignFlowRate';'ZoneMixing';'ZoneControl:Thermostat';'ZoneControl:Humidistat';};
%Known missing & not needed (first few it does implicitly)
not_loaded = {'OutdoorAir:Mixer';'OutdoorAir:Node';'OutdoorAir:NodeList';'Pipe:Adiabatic';...
    'ZoneVentilation:DesignFlowRate';...
    'AirLoopHVAC:OutdoorAirSystem:EquipmentList';...
    'CondenserEquipmentList';'CondenserEquipmentOperationSchemes';...
    'ConnectorList';'ConvergenceLimits';'FuelFactors';'GlobalGeometryRules';...
    'Output:PreprocessorMessage';'OutputControl:Table:Style';...
    'Output:SQLite';'Output:Table:SummaryReports';'OutputControl:ReportingTolerances';...
    'PlantEquipmentList';'PlantEquipmentOperationSchemes';...
    'PlantEquipmentOperation:CoolingLoad';'PlantEquipmentOperation:HeatingLoad';...
    'RunPeriodControl:DaylightSavingTime';'ScheduleTypeLimits';'SimulationControl';...
    'Shading:Building:Detailed';'Shading:Zone:Detailed';'ShadowCalculation';...
    'UtilityCost:Charge:Block';'UtilityCost:Charge:Simple';'UtilityCost:Qualify';...
    'UtilityCost:Tariff';'UtilityCost:Variable';'ZoneAirHeatBalanceAlgorithm'};
%Need an example in idf to include
need_example = {'AirTerminal:SingleDuct:Mixer';'AirLoopHVAC:UnitaryHeatPump:AirToAir';...
    'AvailabilityManager:ScheduledOff';'AvailabilityManager:NightVentilation';...
    'AvailabilityManager:DifferentialThermostat';'AvailabilityManager:HighTemperatureTurnOff';...
    'AvailabilityManager:HighTemperatureTurnOn';'AvailabilityManager:LowTemperatureTurnOff';...
    'AvailabilityManager:LowTemperatureTurnOn';'AvailabilityManager:HybridVentilation';...
    'AvailabilityManager:ScheduledOn';...
    'Chiller:ConstantCOP';'Chiller:Absorption';'Chiller:Absorption:Indirect';...
    'Chiller:CombustionTurbine';'ChillerHeater:Absorption:DirectFired';...
    'ChillerHeater:Absorption:DoubleEffect';'Chiller:EngineDriven';...
    'CoolingTower:TwoSpeed';'CoolingTower:VariableSpeed:Merkel';...
    'Coil:Cooling:WaterToAirHeatPump:VariableSpeedEquationFit';...
    'LoadProfile:Plant';'ZoneHVAC:WindowAirConditioner';'ZoneHVAC:PackagedTerminalHeatPump';};
A = unique(objects.type);
B = A(~ismember(A,loaded) & ~ismember(A,not_loaded));
if ~isempty(B)
    if any(ismember(B,need_example))
        disp('idf has a new equipment type I need to include')
    else
        disp('missing object types during import')
    end
end
end%Ends function import_idf

function [objects,ground_temps,convection] = read_objects(filename)
%reads in the text file into 'objects'
%only exception to pattern is ground temperatures (single line)
if isempty(strfind(filename,'.idf'))
    filename = strcat(filename,'.idf');
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

i = 0;
j = 0;
ground_temps = [];
convection = {};
while i<n_l-1
    i = i+1;
    com = strfind(lines{i},',');
    if ~isempty(com) && isempty(strfind(lines{i},'!'))
        if isempty(lines{i+1})
            if ~isempty(strfind(lines{i},'Site:GroundTemperature:BuildingSurface,'))
                com = strfind(lines{i},',');
                com(end+1) = strfind(lines{i},';');
                ground_temps = zeros(12,1);
                for k = 1:1:12
                    ground_temps(k) = str2double(lines{i}(com(k)+1:com(k+1)-1));
                end
            elseif ~isempty(strfind(lines{i},'SurfaceConvectionAlgorithm:'))
                col = strfind(lines{i},':');
                com = strfind(lines{i},',');
                if strcmpi(lines{i}(col+1:com-1),'Inside')
                    convection.interior = lines{i}(com+1:end-1);
                else
                    convection.exterior = lines{i}(com+1:end-1);
                end
            end
        else
            j = j+1;
            objects.type{j,1} = rmv_spaces(lines{i}(1:com(1)-1));
            com = strfind(lines{i+1},',');
            if isempty(com)
                com = strfind(lines{i+1},';');
            end
            objects.name{j,1} = rmv_spaces(lines{i+1}(1:com(1)-1));
            i = i+2;
            exc = strfind(lines{i},'!');
            value = {};
            parameter = {};
            k = 0;
            while ~isempty(exc)
                k = k+1;
                com = strfind(lines{i},',');
                com = [com,strfind(lines{i},';')];
                com = com(com<exc);
                value{k,1} = rmv_spaces(lines{i}(1:com(end)-1));
                parameter{k,1} = rmv_spaces(lines{i}(exc+2:end));
                i = i+1;
                exc = strfind(lines{i},'!');
            end
            objects.value{j,1} = value;
            objects.parameter{j,1} = parameter;
        end
    end
end
end%Ends function read_objects

function obj = load_object(objects,cat,type,prop_list,n_repeat)
%given objects of a certain type, and the parameter names for which there
%should be values, create a structure with the values/strings from the idf
j1 = 0;
obj = [];
for j = 1:1:length(type)
    if isempty(cat)
        A = f_index(type{j},objects.type);
    else
        if isempty(type{j})
            A = f_index(cat,objects.type);
        elseif strcmp(cat,'ZoneHVAC') && strcmp(type{j},'UnitaryHeatCool')%exception to rule, because I want it lumped in with ZoneHVAC category
            A = f_index('AirLoopHVAC:UnitaryHeatCool',objects.type);%Should seperate this back out.
        else
            A = f_index(strcat(cat,':',type{j}),objects.type);
        end
        obj.type(j1+1:j1+length(A),1) = type(j);
    end
    n = length(A);
    obj.name(j1+1:j1+n,1) = objects.name(A);
    if n>0
        prop = prop_list{min(j,length(prop_list))};
        np = length(prop);
        for i = 1:1:length(prop)
            for k = 1:1:n
                value = objects.value{A(k)};
                if length(value)<i
                    value(i) = {''};
                end
                if ~isfield(obj,prop{i})
                    obj.(prop{i}) = cell(length(obj.name),1);
                    obj.(prop{i})(1:j1) = {''};
                end
                if i>(np - n_repeat)
                    obj.(prop{i})(j1+k,1) = {value(i:n_repeat:end)};
                else
                    if isnumeric(obj.(prop{i}))
                        obj.(prop{i})(j1+k,1) = str2double(value{i});
                    else
                        if ~isempty(strfind(prop{i},'schedule'))
                            obj.(prop{i})(j1+k,1) = {fix_name(value{i})};
                        else
                            obj.(prop{i})(j1+k,1) = value(i);
                        end
                    end
                end
            end   
            if i<=(np - n_repeat)
                %try converting to numeric values from strings
                if ~isnumeric(obj.(prop{i})) && isempty(strfind(prop{i},'schedule'))
                    numbers = str2double(obj.(prop{i}));
                    if any(~isnan(numbers))
                        obj.(prop{i}) = numbers;
                    elseif all(ismember(obj.(prop{i}),{'';'AUTOSIZE';'autocalculate';'AutoSize';'autosize';'Autosize';}))
                        obj.(prop{i}) = nan(length(obj.name),1);
                    end
                end
            end
        end
        if j>1%make sure all inputs have same length when objects of different types have different inputs.
            obj_f = fieldnames(obj);
            obj_f = obj_f(~ismember(obj_f,[prop;'name';'type';]));
            for i = 1:1:length(obj_f)
                if isnumeric(obj.(obj_f{i}))
                    obj.(obj_f{i})(j1+1:j1+n,1) = NaN;
                else
                    obj.(obj_f{i})(j1+1:j1+n,1) = {''};
                end
            end
        end
        j1 = j1+n;
    end
end
end%ends function load_object

function name = fix_name(name)
%matlab doesn't like spaces, commas, dashes, etc in structure names
name = rmv_spaces(name);
if ~isempty(name)
    name = strrep(strrep(name,' ','_'),',','');
    name = strrep(strrep(name,'-','_'),'/','_');
    name = strrep(strrep(name,':','_'),'#','number');
    name = strrep(strrep(name,'(',''),')','');
    if length(name)>60
        name = name(1:60);
    end
end
end%Ends function fix_name

function name = rmv_spaces(name)
%removes spaces at begining and end of line as this can vary based on format of idf
while ~isempty(name) && strcmp(name(1),' ')
    name = name(2:end);
end
while ~isempty(name) && strcmp(name(end),' ')
    name = name(1:end-1);
end
end%Ends function rmv_spaces

function date = import_date(val,ts)
%convert simulation period (month/day) and start day of week into a date vector
%ts, 'timestep' is in hours
switch val{5}
    case 'Monday'
        Y = 2018;
    case 'Sunday'
        Y = 2017;
    case 'Saturday'
        Y = 2011;
    case 'Friday'
        Y = 2010;
    case 'Thursday'
        Y = 2015;
    case 'Wednesday'
        Y = 2014;
    case 'Tuesday'
        Y = 2013;
    otherwise
            disp('need other starting days for simulation')
end
D1 = datenum([Y,str2double(val{1}),str2double(val{2}),0,0,0]);
D2 = datenum([Y,str2double(val{3}),str2double(val{4}),0,0,0])+1;
date = linspace(D1,D2,(D2-D1)*24/ts+1)';
end%Ends function import_date

function occ = import_occupancy(names,values,zone)
for k = 1:1:length(names)
    value = values{k};
    occ.name(k,1) = names(k);
    z = f_index(value{1},zone.name);
    occ.zone(k,1) = z;
    occ.schedule(k,1) = {fix_name(value{2})};
    if strcmpi(value{3},'Area/Person')
        occ.nominal(k,1) = 1/str2double(value{6});
    elseif strcmpi(value{3},'People/Area')
        occ.nominal(k,1) = str2double(value{5});
    elseif strcmpi(value{3},'People')
        people = str2double(value{4});
        occ.nominal(k,1) = people/zone.floor_area(z);
    else
        disp('need more categories in import_occupancy')
    end
    occ.radiant(k,1) = str2double(value{7});
    occ.convected(k,1) = 1 - occ.radiant(k);
    occ.activity(k,1) = {fix_name(value{9})};
    occ.work_eff(k,1) = {fix_name(value{14})};
end
end%Ends function import_occupancy

function light = import_lighting(names,values,zone,occ)
for k = 1:1:length(names)
    value = values{k};
    light.name(k,1) = names(k);
    z = f_index(value{1},zone.name);
    light.zone(k,1) = z;
    light.schedule(k,1) = {fix_name(value{2})};
    if strcmpi(value{3},'LightingLevel')
        lighting = str2double(value{4});
    elseif strcmpi(value{3},'Watts/Area')
        lighting = str2double(value{5})*zone.floor_area(z);
    elseif strcmpi(value{3},'Watts/Person')
        k2 = f_index(light.zone(k),occ.zone);
        lighting = occ.nominal(k2)*str2double(value{6});
    else
        disp('need more categories in mport_lighting')
    end
    light.nominal(k,1) = lighting/zone.floor_area(z);
    light.radiant(k,1) = str2double(value{8});
    light.visible(k,1) = str2double(value{9});
    light.convected(k,1) = 1 - light.radiant(k) - light.visible(k);
end
end%Ends function import_lighting

function equip = import_equipment(names,values,zone,occ)
equip = [];
for k = 1:1:length(names)
    value = values{k};
    equip.name(k,1) = names(k);
    z = f_index(value{1},zone.name);
    equip.zone(k,1) = z;
    equip.schedule(k,1) = {fix_name(value{2})};
    if strcmpi(value{3},'EquipmentLevel')
        power = str2double(value{4});
    elseif strcmpi(value{3},'Watts/Area')
        power = str2double(value{5})*zone.floor_area(z);
    elseif strcmpi(value{3},'Watts/Person')
        k2 = f_index(equip.zone(k),occ.zone);
        power = occ.nominal(k2)*str2double(value{6});
    else
        disp('need more categories in import_quipment')
    end
    equip.nominal(k,1) = power/zone.floor_area(z);
    equip.latent(k,1) = str2double(value{7});
    equip.radiant(k,1) = str2double(value{8});
    equip.lost(k,1) = str2double(value{9});
    equip.convected(k,1) = 1 - equip.latent(k) - equip.radiant(k) - equip.lost(k);
end
end%Ends function import_equipment

function infiltration = import_infiltration(objects,zone)
%outdoor air infiltrating zone through cracks
infiltration.nominal = [];
A = f_index('ZoneInfiltration:DesignFlowRate',objects.type);
for j = 1:1:length(A)
    value = objects.value{A(j)};
    z = f_index(value{1},zone.name);%zone associated with infiltration
    infiltration.zone(j,1) = z;
    infiltration.schedule(j,1) = {fix_name(value{2})};
    switch value{3}
        case 'Flow/Zone'
            infiltration.nominal(j,1) = str2double(value{4});
        case 'Flow/Area'
            infiltration.nominal(j,1) = zone.floor_area(z)*str2double(value{5});
        case 'Flow/ExteriorArea'
            infiltration.nominal(j,1) = zone.exterior_area(z)*str2double(value{6});
        case 'AirChanges/Hour'
            infiltration.nominal(j,1) = zone.volume(z)/3600*str2double(value{7});
    end
end
end%Ends function import_infiltration

function mixing = import_mixing(objects,zone,schedule)
%mixing between zones
mixing.nominal = [];
A = f_index('ZoneMixing',objects.type);
for j = 1:1:length(A)
    value = objects.value{A(j)};
    mixing.receiving_zone(j,1) =  f_index(value{1},zone.name);%zone associated with infiltration
    mixing.schedule(j,1) = {fix_schedule_name(fix_name(value{2}),schedule)};
    mixing.type(j,1) = value(3);
    if strcmpi(value{3},'Flow/Zone')
            mixing.nominal(j,1) = str2double(value{4});
    elseif strcmpi(value{3},'Flow/Area')
            mixing.nominal(j,1) = zone.floor_area(mixing.receiving_zone(j))*str2double(value{5});
    elseif strcmpi(value{3},'Flow/Person')
            mixing.nominal(j,1) = str2double(value{6});
    elseif strcmpi(value{3},'AirChanges/Hour')
            mixing.nominal(j,1) = zone.volume(mixing.receiving_zone(j))/3600*str2double(value{7});
    else
        disp('need more categories in zoneMixing')
    end
    mixing.source_zone(j,1)  = f_index(value{8},zone.name);
end
end%Ends function import_mixing

function [zone,surf,window] = zone_geometry(zone,surf,window)
n_z = length(zone.name);
n_sur = length(surf.name);
n_w = length(window.name);
ceil_area = zeros(n_z,1); 
floor_height = zeros(n_z,1);
floor_area = zeros(n_z,1);
ceil_height = zeros(n_z,1);
for i = 1:1:length(surf.name)
    k = f_index(surf.zone{i},zone.name);
    if strcmp(surf.type{i},'Floor')
        floor_area(k) = floor_area(k) + surf.area(i);
        floor_height(k) = mean(surf.vertices{i}(:,3));%z-coordinate
    elseif any(strcmp(surf.type{i},{'Roof';'Ceiling'}))
        ceil_area(k) = ceil_area(k) + surf.area(i);
    end
end
surf.zone_index = zeros(n_sur,1);
zone.exterior_area = zeros(n_z,1);
zone.windows = cell(n_z,1);
zone.surfaces = cell(n_z,1);
zone.origin(:,1) = zone.x_origin;
zone.origin(:,2) = zone.y_origin;
zone.origin(:,3) = zone.z_origin;
for i = 1:1:n_sur
    k = f_index(surf.zone{i},zone.name);
    surf.zone_index(i) = k;
    zone.surfaces(k) = {[zone.surfaces{k};i]};
    if ~strcmp(surf.type{i},'InternalMass')
        surf.vertices{i} = surf.vertices{i} + ones(length(surf.vertices{i}(:,1)),1)*zone.origin(k,:);%change 3-d position of surface based on zone "origin"
    end
    if any(strcmp(surf.type{i},{'Roof';'Ceiling'}))
        ceil_height(k) = ceil_height(k) + mean(surf.vertices{i}(:,3))*surf.area(i)/ceil_area(k);%z-coordinate 
    end
    if length(surf.boundary)>=i && strcmp(surf.boundary{i},'Outdoors')
        zone.exterior_area(k) = zone.exterior_area(k) + surf.area(i);
    end
end
ceil_height = max(0,ceil_height-floor_height);
volume = ceil_height.*floor_area;
%load in items that were not specified in idf
zone.floor_area(isnan(zone.floor_area)) = floor_area(isnan(zone.floor_area));
zone.ceil_height(isnan(zone.ceil_height)) = ceil_height(isnan(zone.ceil_height));
zone.volume(isnan(zone.volume)) = volume(isnan(zone.volume));

for i = 1:1:length(window.name)
    s = f_index(window.surf_name{i},surf.name);
    window.surface(i,1) = s;
    window.zone(i,1) = surf.zone(s);
    k = f_index(window.zone{i},zone.name);
    zone.windows(k) = {[zone.windows{k};i]};
    window.zone_index(i,1) = k;
    vert = length(window.vertices{i}(:,1));
    window.vertices{i} = window.vertices{i} + ones(vert,1)*zone.origin(k,:);
    surf.area(window.surface(i)) = surf.area(window.surface(i)) - window.area(i);
    window.elevation(i,1) =  mean(window.vertices{i,1}(:,3)) - min(surf.vertices{s}(:,3));%height from floor
    if window.normal(i,3)==1 || strcmp(surf.type{window.surface(i)},'Roof')
        window.elevation(i,1) = window.elevation(i,1) + zone.ceil_height(window.zone_index(i,1));
    end
end

zone.view_factors = zeros(n_sur+n_w);
surf.area_rad = surf.area;%in case some internal mass objects surface area must be 'reduced' to enforce reciprocity (area of internal mass greater than area of surfaces in zone)
im = f_index('InternalMass',surf.type);
surf.area_rad(im) = surf.area_rad(im)/2;
for k = 1:1:n_z
    %compute internal view factors in single zone
    sur_k = f_index(k,surf.zone_index);
    win_k = f_index(k,window.zone_index);
    im_sur = f_index('InternalMass',surf.type(sur_k));%
    vf_zone = calc_view_factors(surf.area_rad(sur_k),window.area(win_k),im_sur,[surf.normal(sur_k,:);window.normal(win_k,:)]);
    %%combine all view factors into 1 matrix [surfaces;windows];
    sur = [sur_k;win_k+n_sur];
    for j = 1:1:length(sur)
        zone.view_factors(sur(j),sur) = vf_zone(j,:);
        zone.view_factors(sur,sur(j)) = vf_zone(:,j);
    end
end

ext = sort([f_index('Outdoors',surf.boundary);f_index('Ground',surf.boundary)]);
surf.exterior.name = surf.name(ext,:);
surf.exterior.normal = surf.normal(ext,:);
surf.exterior.cos_phi = surf.cos_phi(ext);
surf.exterior.thermal_absorptance = surf.absorptance.thermal(ext,2);
surf.exterior.solar_absorptance = surf.absorptance.solar(ext,2);
surf.exterior.boundary = surf.boundary(ext);
surf.exterior.roughness_factor = surf.roughness_factor(ext);
surf.exterior.perimeter = surf.perimeter(ext);
surf.exterior.area = surf.area(ext);%area without windows
end%Ends function zone_geometry

function [vf,area] = calc_view_factors(area_s,area_w,im_sur,s_norm)
n = length(area_s);
w = length(area_w);
vf = zeros(n+w,n+w);
area = [area_s;area_w];
net_area = sum(area);
for i = 1:1:n+w
    seen_area = net_area;% - area(i);
    if any(i == im_sur)
        seen_area = seen_area - area(i);
    else
        for j = 1:1:n+w
            if ~any(j == im_sur) && dot(s_norm(i,:),s_norm(j,:))>.985 %parallel surfaces (within 10 degrees) = zero view factor, pg 103 of manual
                seen_area = seen_area - area(j);
            end
        end
    end
    for j = 1:1:n+w
        if i~=j && (dot(s_norm(i,:),s_norm(j,:))<=.985 || any(i == im_sur) || any(j == im_sur))
            vf(i,j) = area(j)/seen_area;
            %otherwise do nothing window/surfaces dont see each other
        end
    end
end
%check reciprocity!: 
area_i = area*ones(1,n+w);
vf_error = ones(n+w,n+w);
while any(any(abs(vf_error)>1e-6))
    vf_scaled = vf.*area_i;
    vf_error = (vf_scaled - vf_scaled');
    perc_error = vf_error./max(vf_scaled,vf_scaled');
    perc_error(isnan(perc_error)) = 0;
    perc_error(isinf(perc_error)) = 0;
    vf = vf.*(1-.5*perc_error);
    vf = vf.*((1./sum(vf,2))*ones(1,n+w));%ensure all add up to 1
end
end%Ends function calc_view_factors

function outdoor_air = import_idf_outdoor_air(names,values)
for k = 1:1:length(names)
    outdoor_air.name(k,1) = names(k);
    outdoor_air.method(k,1) = values{k}(1);
    switch values{k}{1}
        case 'Flow/Person'
            j = 2;
        case 'Flow/Area'
            j = 3;
        case 'Flow/Zone'
            j = 4;
    end
    val = str2double(values{k}{j});
    if isnan(val)
        outdoor_air.value(k,1) = 0;
    else
        outdoor_air.value(k,1) = val;
    end
end
end%Ends function import_idf_outdoor_air

function surf = import_idf_surfaces(objects,construct)
A = f_index('BuildingSurface:Detailed',objects.type);
A2 = f_index('InternalMass',objects.type);
prop1 = {'type';'construct';'zone';'boundary';'object';};%'sun';'wind';'view2ground';'vertices'};
for j = 1:1:length(A)
    value = objects.value{A(j)};
    surf.name(j,1) = objects.name(A(j));
    for p = 1:1:length(prop1)
        surf.(prop1{p})(j,1) = value(p);
    end
    if strcmp(surf.boundary(j),surf.name(j))
        surf.adiabatic(j,1) = true;  
    else
        surf.adiabatic(j,1) = false;  
    end
    surf.vertices(j,1) = {surface_vertices(value(10:9+str2double(value{9})))};
    surf.height(j,1)  = max(surf.vertices{j}(:,3)) - min(surf.vertices{j}(:,3));
    [surf.area(j,1),surf.perimeter(j,1),surf.normal(j,:)] = find_area(surf.vertices{j},surf.type{j});
    c = f_index(surf.construct{j},construct.name);
    surf.roughness(j,1) = construct.roughness(c);%outside surface material roughness
    surf.absorptance.solar(j,:) = construct.solar(c,:);%absorptance
    surf.absorptance.thermal(j,:) = construct.thermal(c,:);%absorptance
    surf.absorptance.visible(j,:) = construct.visible(c,:);%absorptance
    surf.capacitance(j,1) = construct.capacitance(c);
    surf.thermal_conductivity(j,1) = construct.thermal_conductivity(c);
    surf.node_length(j,1) = construct.node_length(c);
    surf.resistance(j,1) = construct.resistance(c);
    surf.i_nodes(j,1) = length(surf.capacitance);%interior wall states per layer minimum is 3, both surfaces and center
end
surf.roughness_factor(strcmpi(surf.roughness,'VeryRough'),1) = 2.17;
surf.roughness_factor(strcmpi(surf.roughness,'Rough'),1) = 1.67;
surf.roughness_factor(strcmpi(surf.roughness,'MediumRough'),1) = 1.52;
surf.roughness_factor(strcmpi(surf.roughness,'MediumSmooth'),1) = 1.13;
surf.roughness_factor(strcmpi(surf.roughness,'Smooth'),1) = 1.11;
for k = 1:1:length(A2)
    value = objects.value{A2(k)};
    surf.name(j+k,1) = objects.name(A2(k));
    surf.type(j+k,1) = {'InternalMass'};
    surf.construct(j+k,1) = value(1);
    surf.zone(j+k,1) = value(2);
    surf.area(j+k,1) = str2double(value{3});
    surf.normal(j+k,:) = [0 0 1]; %facing upward
    surf.vertices(j+k) = {''};
    surf.boundary(j+k,1) = {'Surface'};
    surf.object(j+k,1) = {objects.name(A2(k))};  
    surf.adiabatic(j+k,1) = true;  
    c = f_index(surf.construct{j+k},construct.name);
    surf.absorptance.solar(j+k,:) = construct.solar(c,:);%absorptance
    surf.absorptance.thermal(j+k,:) = construct.thermal(c,:);%absorptance
    surf.absorptance.visible(j+k,:) = construct.visible(c,:);%absorptance
    surf.capacitance(j+k,1) = construct.capacitance(c);
%     surf.capacitance(j+k,1) = {construct.capacitance{c}/4};
    surf.thermal_conductivity(j+k,1) = construct.thermal_conductivity(c);
    surf.node_length(j+k,1) = construct.node_length(c);
    surf.i_nodes(j,1) = length(surf.capacitance);%interior wall states per layer minimum is 3, both surfaces and center
end
surf.cos_phi = surf.normal(:,3)./(surf.normal(:,1).^2 + surf.normal(:,2).^2 + surf.normal(:,3).^2).^.5; %portion of surface pointed upwards

surf.orientation.vertical = abs(surf.normal(:,3))<1e-3; 
surf.orientation.horizontal =  abs(surf.normal(:,1))<1e-3 &  abs(surf.normal(:,2))<1e-3;
surf.orientation.tilted = ~surf.orientation.vertical & ~surf.orientation.horizontal;
surf.orientation.ceiling = surf.orientation.horizontal & surf.normal(:,3)<0;
surf.orientation.floor = surf.orientation.horizontal & surf.normal(:,3)>0;
end%ends function import_idf_surfaces

function window = import_window(objects,construct)
%%Windows & doors
B = f_index('FenestrationSurface:Detailed',objects.type);
w = 0;
window.name = {};
for j = 1:1:length(B)
    value = objects.value{B(j)};
    c = f_index(value{2},construct.name);
    if strcmpi(value(1),'Window') %length(c) == 1
        w = w+1;
        window.name(w,1) = objects.name(B(j));
        window.type(w,1) = construct.type(c);
        window.construct(w,1) = value(2);
        window.surf_name(w,1) = value(3);
        window.frame_name(w,1) = value(7);
        window.vertices(w,1) = {surface_vertices(value(10:9+str2double(value{9})))};
        [window.area(w,1),window.perimeter(w,1),window.normal(w,:),window.height(w,1)] = find_area(window.vertices{w,1},'');
        prop = {'normal_transmittance';'normal_reflectance';'emittance';'thermal_conductivity';'solar_heat_gain';'gap_thickness';'transmittance';'reflectance';'thickness';};
        for p = 1:1:length(prop)
            window.(prop{p})(w,:) = construct.(prop{p})(c,:);
        end
    end
end
window.cos_phi = window.normal(:,3)./(window.normal(:,1).^2 + window.normal(:,2).^2 + window.normal(:,3).^2).^.5; %portion of surface pointed upwards
window.orientation.vertical = abs(window.normal(:,3))<1e-3; 
window.orientation.horizontal =  abs(window.normal(:,1))<1e-3 &  abs(window.normal(:,2))<1e-3;
window.orientation.tilted = ~window.orientation.vertical & ~window.orientation.horizontal;
window.orientation.ceiling = window.orientation.horizontal & window.normal(:,3)<0;
window.orientation.floor = window.orientation.horizontal & window.normal(:,3)>0;
end%Ends function import_window

function door = import_door(objects,construct)
B = f_index('FenestrationSurface:Detailed',objects.type);
d = 0;
door.name = {};
for j = 1:1:length(B)
    value = objects.value{B(j)};
    c = f_index(value{2},construct.name);
    if length(c) == 1
        d = d+1;
        door.name(d,1) = objects.name(B(j));
        disp('need to be able to import doors and replace portion of wall surface')
    end
end
end%Ends function import_door

function vertices = surface_vertices(vert)
vertices = zeros(length(vert),3);
for j = 1:1:length(vert)
    com = strfind(vert{j},',');
    vertices(j,1) = str2double(rmv_spaces(vert{j}(1:com(1)-1)));
    vertices(j,2) = str2double(rmv_spaces(vert{j}(com(1)+1:com(2)-1)));
    vertices(j,3) = str2double(rmv_spaces(vert{j}(com(2)+1:end)));
end
end%Ends function surface_vertices

function [construct,window_construct] = read_construction(names,value,materials,window_materials,i_nodes)
s = 0;
w = 0;
for k = 1:1:length(value)
    mat = [];
    for i = 1:1:length(value{k})
        m = f_index(value{k}{i},materials.name);
        if ~isempty(m)
            mat(end+1) = m;
        end
    end
    if length(mat) == length(value{k})
        s = s+1;
        construct.name(s,1) = names(k);
        construct.roughness(s,1) = materials.roughness(mat(1));%outside surface material roughness
        [construct.thermal_conductivity(s,1),construct.node_length(s,1),construct.capacitance(s,1),construct.resistance(s,1)] = material_conductivity(mat,materials,i_nodes);
        [construct.solar(s,:),construct.thermal(s,:),construct.visible(s,:)] = absorptance(mat,materials);
    else %window
        w = w+1;
        for i = 1:1:length(value{k})
            m = f_index(value{k}{i},window_materials.name);
            if ~isempty(m)
                mat(end+1) = m;
            end
        end
        window_construct.name(w,1) = names(k);
        window_construct.type(w,1) = window_materials.type(mat(1));
        window_construct = window_glazing(window_construct,window_materials,mat,w);
    end
end
end%ends function import_idf_construction

function racks = rack_cases(racks,cases,names,values)
racks.cases = cell(length(racks.name),1);
for i = 1:1:length(names)
    r = f_index(names{i},racks.case_name);
    value = values{i};
    c_num = [];
    for k = 1:1:length(value)
        c_num(end+1) = f_index(value{k},cases.name);
    end
    racks.cases(r) = {c_num};
end
for i = 1:1:length(racks.cases)
    if isempty(racks.cases{i})
        racks.cases(i) = {f_index(racks.case_name{i},cases.name)};
    end
end
end%Ends function rack_cases

function sched = import_schedules(names,value)
for k = 1:1:length(value)
    name = fix_name(names{k});
    sched.(name).type = value{k}{1};
    sched.(name).seasons = [];
    sched.(name).ramp = 1e-4;
    i = 2;
    n_v = length(value{k});
    season = 0;
    while i<=n_v
        if strcmpi(value{k}{i}(1:7),'Through')
            bksls = strfind(value{k}{i},'/');
            month = str2double(rmv_spaces(value{k}{i}(bksls(1)-2:bksls(1)-1)));
            day = str2double(rmv_spaces(value{k}{i}(bksls(1)+1:end)));
            dy = [31 28 31 30 31 30 31 31 30 31 30 31];
            sched.(name).seasons(end+1) = sum(dy(1:(month-1)))*24 + 24*day;
            season = season + 1;
            sched.(name).interpolate(season) = false;
            i = i+1;
        elseif strcmpi(value{k}{i}(1:3),'For')
            day = rmv_spaces(value{k}{i}(5:end));
            day = break_day(day);
            i = i+1;
            if strcmpi(value{k}{i},'Interpolate:No')
                i = i+1;
            elseif strcmpi(value{k}{i},'Interpolate:Yes')
                sched.(name).interpolate(season) = true;
                i = i+1;
            end
            mat = [];
            j = 0;
            while i<=n_v && strcmpi(value{k}{i}(1:5),'Until')
                col = strfind(value{k}{i},':');
                com = strfind(value{k}{i},',');
                mat(j+1,1) = str2double(value{k}{i}(col(1)+1:col(2)-1)) + str2double(value{k}{i}(col(2)+1:com(1)-1))/60;
                mat(j+1,2) = str2double(value{k}{i}(com(1)+1:end));
                i = i+1;
                j = j+1;
            end
            for j = 1:1:length(day)
                sched.(name).(day{j}){season} = [0,mat(1,2);mat;];
                for s = 1:1:season-1
                    if isempty(sched.(name).(day{j}){s})
                        sched.(name).(day{j}){s} = sched.(name).(day{j}){season};
                    end
                end
            end
        else
            i = i+1;%skip miscelaneous line
        end
    end
end
end%ends function import_schedules

function [K,dX,C,R] = material_conductivity(layers,material,i_nodes)
n = length(layers);
index = (1-i_nodes):0;
K = [];
dX = [];
R = [];
C_total = 0;
for i = 1:1:n
    if isfield(material,'thickness') && ~isnan(material.thickness(layers(i)))
        C_total = C_total +  material.density(layers(i))*material.specheat(layers(i))*material.thickness(layers(i));%capacitance in J/K*m^2
    end
end
C = 0;
for i = 1:1:n
    if isfield(material,'thickness') && ~isnan(material.thickness(layers(i)))
        
        C_l = material.density(layers(i))*material.specheat(layers(i))*material.thickness(layers(i));%capacitance in J/K*m^2
        if C_l/C_total>0.04
            index = index + i_nodes;
            C(end) = C(end) + .5*C_l/i_nodes;
            C(index+1) = C_l/i_nodes;
            C(end) = C(end) - .5*C_l/i_nodes;
            if length(R) == index(1)
                R1 = R(index(1));
            else R1 = 0;
            end
            dX(index) = material.thickness(layers(i))/i_nodes;%m
            R(index) = material.thickness(layers(i))/material.conductivity(layers(i))/i_nodes;% K*m^2/W 
            R(index(1)) = R1 + R(index(1));
            K(index) = dX(index)./R(index);% W/m*K
        else%create a single node
            index = index + 1;
            C(end:end+1) = [C(end),0] + .5*C_l;
            R(end+1) = material.thickness(layers(i))/material.conductivity(layers(i));% K*m^2/W 
            dX(end+1) = material.thickness(layers(i));
            K(end+1) = material.conductivity(layers(i));% W/m*K   
%             %don't create more nodes, lump in with portion of surface with greater capacitance
%             C(end) = C(end) + C_l;
%             R(end) = R(end) + material.thickness(layers(i))/material.conductivity(layers(i));% K*m^2/W 
%             K(end) = dX(end)/R(end);% W/m*K            
        end
    elseif isfield(material,'thermal_resistance') && ~isnan(material.thermal_resistance(layers(i)))% no mass layers
        if C_total == 0
            %air wall, find thickness that matches R and use that for capacitance
            k_air = 2.5e-2; %W/m*K
            rho_air = 1.204;%kg/m^3
            Cp_air = 1006;%Specific Heat {J/kg-K}
            t = material.thermal_resistance(layers(i))*k_air;%thickness in m
            air_cap = t*rho_air*Cp_air; %capacitance in J/K*m^2
            K(i) = 0.0262;%air equivelant, air = 0.0262 W/m*K
            dX(i) = material.thermal_resistance(layers(i))./K(i);%air equivelant, air = 0.0262 W/m*K
            R(i) =  material.thermal_resistance(layers(i));% K*m^2/W 
            C(i:i+1,1) = air_cap/2;
        else%no-mass layer lumped into previous node
            if i == 1
                R(1) = material.thermal_resistance(layers(i));% K*m^2/W       
            else
                R(end) = R(end) + material.thermal_resistance(layers(i));% K*m^2/W 
                K(end) = dX(end)/R(end);% W/m*K      
            end
        end
        
    end
end
% %average
% C = sum(C)*[.125,.25,.25,.25,.125];
% R = sum(1./(K./dX));%resistances in series add up
% dX = sum(dX)*[.25,.25,.25,.25];
% K = (1./(R*[.25,.25,.25,.25])).*dX;
% R = sum(R)*[.25,.25,.25,.25];%resistances in series add up
K = {K};
dX = {dX};
C = {C};
R = {R};
end%Ends function material_conductivity

function [solar,thermal,visible] = absorptance(layers,material)
db = ones(1,2);
solar = db*material.solar_absorptance(layers(end));%inside surface material;
thermal = db*material.thermal_absorptance(layers(end));%inside surface material;
visible = db*material.visible_absorptance(layers(end));%inside surface material;
if length(layers)>1%record properties of outside layer if different
    solar(1,2) = material.solar_absorptance(layers(1));%outside surface material
    thermal(1,2) = material.thermal_absorptance(layers(1));%outside surface material
    visible(1,2) = material.visible_absorptance(layers(1));%outside surface material
end
end%Ends function absorptance

function [area,perimeter,normal,max_height] = find_area(vertices,type)
P1 = vertices(1,:);
P2 = vertices(2,:);
P3 = vertices(3,:);
A = (sum((P2-P1).^2))^.5; %length from P1 to P2
B = (sum((P3-P2).^2))^.5; %length from P2 to P3
C = (sum((P1-P3).^2))^.5; %length from P3 to P1
j = 3;
while (A^2+B^2 +C^2) == (A+B+C)^2 %change point 3 if vertices are on same line
    j = j+1;
    P3 = vertices(j,:);
    B = (sum((P3-P2).^2))^.5; %length from P2 to P3
    C = (sum((P1-P3).^2))^.5; %length from P3 to P1
end
normal_v = [0,0,0];
for j = 1:1:length(vertices(:,1))-2
    normal_v = normal_v + cross((vertices(j+1,:)-vertices(j,:)),(vertices(j+2,:)-vertices(j+1,:)));
end
normal_v2 = cross((P2-P1),(P3-P2));
a = normal_v(1);%x-value of normal vector 
b = normal_v(2);%y-value of normal vector 
c = normal_v(3);%z-value of normal vector 
mag = (a^2 + b^2 + c^2)^.5;%magnitude of normal vector
normal = normal_v/mag;
normal2 = normal_v2/norm(normal_v2);
if any(abs(normal - normal2)>1e-4)%trouble getting correct normal vector direction with non-convex floor or ceiling shape
    if strcmp(type,'Floor')
        normal = [0,0,-1];
    else
        normal = [0,0,1];
    end
end
n_v = length(vertices(:,1));
if a~=0 || b~=0 %not in X-Y plane, need to rotate
    if dot([a,b,0],[1,0,0]) == 0 %single rotation about x-axis
        theta_z = sign(normal(2))*acos(dot(normal,[0,0,1]));%angle between normal and z-axis (rotate around y-axis)
        R_x = [1 0 0; 0 cos(theta_z) -sin(theta_z); 0 sin(theta_z) cos(theta_z);];
        for k = 1:1:n_v
            vertices(k,:) = (R_x*vertices(k,:)')';
        end
    else
        theta_x = -sign(normal(2))*acos(dot([a,b,0],[1,0,0])/(a^2+b^2)^.5);%angle between shadow and x-axis (rotate around z-axis)
        R_z = [cos(theta_x) -sin(theta_x) 0; sin(theta_x) cos(theta_x) 0; 0 0 1];
        normal2 = (R_z*normal')';
        theta_z = -sign(normal(1))*acos(dot(normal2,[0,0,1]));%angle between normal and z-axis (rotate around y-axis)
        R_y = [cos(theta_z) 0 sin(theta_z); 0 1 0; -sin(theta_z) 0 cos(theta_z);];
        for k = 1:1:n_v
            vertices(k,:) = (R_y*R_z*vertices(k,:)')';
        end
    end
end
area = abs(sum((vertices([2:n_v 1],1) - vertices(:,1)).* (vertices([2:n_v 1],2) + vertices(:,2)))/2);%area under each segment in x-dir, length of segment * average height in y
perimeter = abs(sum(((vertices([2:n_v 1],1) - vertices(:,1)).^2 + (vertices([2:n_v 1],2) - vertices(:,2)).^2).^.5));
dist = zeros(n_v,n_v);
for i = 1:1:n_v
    for j = i+1:1:n_v
        dist(i,j) = ((vertices(j,1) - vertices(i,1))^2 + (vertices(j,2) - vertices(i,2))^2 + (vertices(j,3) - vertices(i,3))^2)^.5;
    end
end
max_height = max(max(dist));
end%Ends function find_area

function name = fix_schedule_name(name,schedule)
f = fieldnames(schedule);
if ~any(strcmp(name,f))
    k = f_index(name,f);
    if ~isempty(k)
        name = f{k};
    end
end
end%Ends function fix_schedule_name

function day = break_day(day)
spac = strfind(day,' ');
if ~isempty(spac)
    all = day;
    day = cell(length(spac)+1,1);
    day(1) = {all(1:spac(1)-1)};
    for i = 1:1:length(spac)-1
        day(i+1) = {all(spac(i)+1:spac(i+1)-1)};
    end
    day(end) = {all(spac(end)+1:end)};
else
    day = {day};
end
end%ends function break_day

function [loop,components] = import_plant(objects)
%% Plant Loops 
prop = {'fluid';'user_fluid';'operation_scheme';'setpoint_node';'max_temperature';'min_temperature';'max_flow_rate';'min_flow_rate';'loop_volume';'supply_inlet';'supply_outlet';'supply_branch_list';'supply_connector_list';'demand_inlet';'demand_outlet';'demand_branch_list';'demand_connector_list';'load_scheme';};
loop = load_object(objects,'',{'PlantLoop';'CondenserLoop'},{prop,prop},0);
loop.loop_order = (1:length(loop.name))';

if ~isempty(loop.loop_order)
    for i = 1:1:length(loop.name)
        if strcmpi(loop.fluid{i},'Water')
            loop.fluid_density(i,1) = 997;%kg/m^3
        else
            disp('need density for plant loop fluid other than water')
        end
    end
    B = f_index('Sizing:Plant',objects.type);
    for k = 1:1:length(B)
        i = f_index(objects.name{B(k)},loop.name);
        if ~isempty(i)
            value = objects.value{B(k)};
            loop.type(i,1) = value(1);
            loop.exit_temperature(i,1) = str2double(value{2});
            loop.temperature_difference(i,1) = str2double(value{3});
        end
    end

    branch_list = load_object(objects,'',{'BranchList';},{{'branches'}},1);
    for i = 1:1:length(loop.name)
        bl_d = f_index(loop.demand_branch_list{i},branch_list.name);
        bl_s = f_index(loop.supply_branch_list{i},branch_list.name);
        loop.branches(i,1) = {[branch_list.branches{bl_d};branch_list.branches{bl_s}]};
    end

    branches = load_object(objects,'',{'Branch';},{{'pressure_curve';'component_type';'component_name';'component_inlet';'component_outlet'}},4);
    branches.loop = zeros(length(branches.name),1);
    for j = 1:1:length(branches.name)
        for l = 1:1:length(loop.branches)
            l_b = f_index(branches.name{j},loop.branches{l});
            if ~isempty(l_b)
                branches.loop(j) = l;
                break
            end
        end
        branches.inlet(j,1) = branches.component_inlet{j}(1);
        branches.outlet(j,1) = branches.component_outlet{j}(end);
    end
    
    splitters = load_object(objects,'Connector',{'Splitter';},{{'inlet';'outlets'}},1);
    splitters = branch_2_node(splitters,branches,'inlet','outlet');
    mixers = load_object(objects,'Connector',{'Mixer';},{{'outlet';'inlets'}},1);
    mixers = branch_2_node(mixers,branches,'outlet','inlet');

    components = read_components(branches,mixers,splitters);
else
    components.name = {};
    loop.name = {};
end
end%Ends function import_plant

function connect = branch_2_node(connect,branches,a,b)
c = strcat(b,'s');
for i = 1:1:length(connect.name)
    b_num = f_index(connect.(a){i},branches.name);
    connect.(a)(i) = branches.(b)(b_num);
    ports = connect.(c){i};
    for j = 1:1:length(ports)
        b2 = f_index(ports{j},branches.name);
        ports(j) = branches.(a)(b2);
    end
    connect.(c){i} = ports;
end
end%Ends function branch_2_node

function [nodes,equip] = setup_nodes(components,loop,d_or_s)
%% create a list of nodes in the order they are connected
%% Identify which branch/loop the node is on, and which # component on that branch it is
nodes.name = {};
equip.name = {};
n = 0;
n2 = 0;
component_loaded = false(length(components.name),1);
for pl = 1:1:length(loop.name)
    if strcmp(d_or_s,'d')
        eob = loop.demand_inlet(pl);%eob is end-of-branch
        outlet = loop.demand_outlet{pl};
    else
        eob = loop.supply_inlet(pl);%eob is end-of-branch
        outlet = loop.supply_outlet{pl};
    end
    nodes.inlet_node(pl) = n+1;
    while ~isempty(eob) 
        for i = 1:1:length(eob)%identify if any branch starts at the current end-of-line
            c = f_index(eob{i},components.inlet);
            if ~isempty(c)
                if any(strcmp('Mixer',components.type(c))) && i~=length(eob)
                    %move on to other branches and come back to mixer
                    %can only be one mixer on a loop, so not a problem this way
                    %this ensures other equipment upstream of mixer is simulated first
                else
                    eob(i) = {''};
                    if any(component_loaded(c)==false)
                        break
                    end
                end
            end
        end
        b = unique(components.branch_number(c));
        b_components = f_index(b,components.branch_number);
        for i = 1:1:length(b_components)
            prev_node = f_index(components.inlet(b_components(i)),nodes.name);
            if isempty(prev_node)
                n = n+1;
                nodes.name(n,1) = components.inlet(b_components(i));
                nodes.branch(n,1) = b;
                nodes.loop(n,1) = pl;
                nodes.bypass(n,1) = false;
                prev_node = n;
            end
            prev_equip = f_index(components.name(b_components(i)),equip.name);
            if isempty(prev_equip)
                n2 = n2+1;
                equip.name(n2,1) = components.name(b_components(i));
                equip.type(n2,1) = components.type(b_components(i));
                equip.component_number(n2,1) = b_components(i);
                equip.inlet(n2,1) = {prev_node};
                equip.outlet_name(n2,1) = components.outlet(b_components(i)); 
                equip.loop(n2,1) = pl;
                equip.bypass(n2,1) = false;
                prev_equip = n2;
            end
            if components.bypass(b_components(i))
                nodes.bypass(prev_node) = true;
                equip.bypass(prev_equip) = true;
            end
        end
        if any(strcmp('Mixer',components.type(b_components)))
            m = f_index('Mixer',components.type(b_components));
            e = f_index(components.name{b_components(m(1))},equip.name);
            c_inlets = f_index(equip.name{e},components.name);
            in_nodes = [];
            for j = 1:1:length(c_inlets)
                k = f_index(components.inlet{c_inlets(j)},nodes.name);
                in_nodes(1,end+1) = k;
            end
            equip.inlet(e) = {in_nodes};
        end
        component_loaded(b_components) = true;
        if any(strcmp('Splitter',components.type(b_components)))
            s = f_index('Splitter',components.type(b_components));
            eob(end+1:end+length(s)) = components.outlet(b_components(s));
        else
            eob(end+1) = components.outlet(b_components(end));
        end

        eob = eob(~strcmpi('',eob));
        eob = unique(eob);
        eob = eob(~strcmpi(outlet,eob));
        if isempty(eob) && ~any(strcmpi(outlet,nodes.name))
            n = n+1;
            nodes.name(n,1) = {outlet};
            nodes.branch(n,1) = b;
            nodes.loop(n,1) = pl;
            nodes.outlet_node(pl) = n;
        end
    end
end
equip.outlet = cell(n2,1);
for i = 1:1:length(nodes.name)
    upstream = f_index(nodes.name(i),equip.outlet_name);
    for j = 1:1:length(upstream)
        if strcmp(equip.type(upstream(j)),'Splitter')
            s = f_index(equip.name{upstream(j)},components.name);
            out_nodes = [];
            outlets = components.outlet(s);
            for l = 1:1:length(outlets)
                out_nodes(end+1) = f_index(outlets{l},nodes.name);
            end
            equip.outlet(upstream(j),1) = {out_nodes};
        else
            equip.outlet(upstream,1) = {i};
        end
    end
end
end%ends function setup_nodes

function [hvac,e_list] = import_hvac(objects,zone_names,schedules)
prop = {'controller';'manager';'supply_flow';'branch_list';'connector_list';'supply_inlet';'demand_outlet';'demand_inlet';'supply_outlet';};
hvac.loop = load_object(objects,'',{'AirLoopHVAC'},{prop},0);
branch_list = load_object(objects,'',{'BranchList';},{{'branches'}},1);
for i = 1:1:length(hvac.loop.name)
    bl = f_index(hvac.loop.branch_list{i},branch_list.name);
    hvac.loop.branches(i,1) = branch_list.branches(bl);
end
hvac.design = load_design_size(objects,hvac.loop);
hvac.branches = load_object(objects,'',{'Branch';},{{'pressure_curve';'component_type';'component_name';'component_inlet';'component_outlet'}},4);
hvac.branches.loop = zeros(length(hvac.branches.name),1);
for j = 1:1:length(hvac.branches.name)
    for l = 1:1:length(hvac.loop.branches)
        l_b = f_index(hvac.branches.name{j},hvac.loop.branches{l});
        if ~isempty(l_b)
            hvac.branches.loop(j) = l;
            break
        end
    end
    hvac.branches.inlet(j,1) = hvac.branches.component_inlet{j}(1);
    hvac.branches.outlet(j,1) = hvac.branches.component_outlet{j}(end);
end
hvac.supply_path = load_object(objects,'AirLoopHVAC',{'SupplyPath';},{{'inlet';'component_type';'component_name'}},2);
for k = 1:1:length(hvac.supply_path.name)
    hvac.supply_path.loop(k,1) = f_index(hvac.supply_path.inlet{k},hvac.loop.demand_inlet);
    hvac.branches.inlet(end+1) = {''};% hvac.supply_path.inlet(k);
    hvac.branches.outlet(end+1) = hvac.supply_path.inlet(k);
    hvac.branches.name(end+1) = hvac.supply_path.name(k);
    hvac.branches.loop(end+1) = 0;%hvac.supply_path.loop(k);
end
hvac.return_path = load_object(objects,'AirLoopHVAC',{'ReturnPath';},{{'outlet';'component_type';'component_name'}},2);
for k = 1:1:length(hvac.return_path.name)
    hvac.return_path.loop(k,1) = f_index(hvac.return_path.outlet{k},hvac.loop.demand_outlet);
    hvac.branches.inlet(end+1) = hvac.return_path.outlet(k);
    hvac.branches.outlet(end+1) = {''};% hvac.return_path.outlet(k);
    hvac.branches.name(end+1) = hvac.return_path.name(k);
    hvac.branches.loop(end+1) = 0;%hvac.return_path.loop(k);
end
%Items with varying length (# of connections etc)
hvac.splitters = load_object(objects,'AirLoopHVAC',{'ZoneSplitter';},{{'inlet';'outlets'}},1);
hvac.mixers = load_object(objects,'AirLoopHVAC',{'ZoneMixer';},{{'outlet';'inlets'}},1);
hvac.plenum = load_object(objects,'AirLoopHVAC',{'ReturnPlenum';},{{'zone';'zone_node';'outlet';'induced_air';'inlets'}},1);

prop = {'hvac_equipment';'inlet';'exhaust';'node_name';'return';};
zone_hvac = load_object(objects,'ZoneHVAC',{'EquipmentConnections'},{prop},0);
for p = 1:1:length(prop)
    if isnumeric(zone_hvac.(prop{p}))
        zone_hvac.(prop{p}) = cell(length(zone_hvac.name),1);
    end
end

A = f_index('ZoneHVAC:EquipmentList',objects.type);
k = 0;
e_list.name = {};
for i = 1:1:length(A)
    value = objects.value{A(i)};
    name =  objects.name(A(i));
    ec_list = f_index(name,zone_hvac.hvac_equipment);
    air_return = zone_hvac.return(ec_list);
    for j = 1:1:length(hvac.plenum.name)
        p_list = f_index(air_return,hvac.plenum.inlets{j});
        if ~isempty(p_list)
            air_return = hvac.plenum.outlet(j);%connecting terminal outlet to air node seen by mixer on return path
            break
        end
    end
    z = f_index(zone_hvac.name{ec_list},zone_names);
    for j = 1:4:length(value)
        k = k+1;
        e_list.name(k,1) = name;
        e_list.type(k,1) = value(j);
        e_list.zone(k,1) = z;
        e_list.e_name(k,1) = value(j+1);
        e_list.cool_sequence(k,1) = str2double(value{j+2});
        e_list.heat_sequence(k,1) = str2double(value{j+3});
        e_list.inlet(k,1) = zone_hvac.inlet(ec_list);
        e_list.outlet(k,1) = zone_hvac.exhaust(ec_list);
        e_list.return(k,1) = air_return;
    end
end
%Terminals
hvac.terminals = load_terminals(objects,hvac.splitters,hvac.loop,e_list);
%All HVAc components
hvac.components = read_components(hvac.branches,hvac.mixers,hvac.splitters);
%%add terminals/distribution units/zones to connect demand equipment
for i = 1:1:length(hvac.terminals.name)
    hvac.components.name(end+1,1) = hvac.terminals.name(i);
    hvac.components.type(end+1,1) = {strcat('AirTerminal:SingleDuct:',hvac.terminals.type{i})};
    hvac.components.inlet(end+1,1) = hvac.terminals.inlet(i);
    hvac.components.outlet(end+1,1) = hvac.terminals.return(i);
    hvac.components.bypass(end+1,1) = false;
    hvac.components.branch(end+1,1) = zone_names(hvac.terminals.zone(i));
    hvac.components.branch_number(end+1,1) = length(hvac.branches.name)+hvac.terminals.zone(i);
    hvac.components.loop(end+1,1) = hvac.terminals.loop(i);
end        
%%Outdoor air
prop1 = {'controller_list';'equipment_list';'manager_list';};
oa_sys = load_object(objects,'AirLoopHVAC',{'OutdoorAirSystem'},{prop1},0);
prop1 = {'controller_type';'controller_name';};
oa_list = load_object(objects,'AirLoopHVAC',{'ControllerList'},{prop1},0);
prop1 = {'relief_outlet';'return_air';'mixed_air';'actuator_node';'min_flow';'max_flow';'economizer_control';'economizer_action';'economizer_max_T';'economizer_max_h';'economizer_max_dp';'enthalpy_limit_curve';'economizer_min_T';'lockout';'min_limit_type';'min_air_schedule';'min_frac_schedule';'max_frac_schedule';'mechanical_controller_name';};
hvac.outdoor_air = load_object(objects,'Controller',{'OutdoorAir'},{prop1},0);
for i = 1:1:length(hvac.outdoor_air.name)
    ln = f_index(hvac.outdoor_air.name{i},oa_list.controller_name);
    sys = f_index(oa_list.name{ln},oa_sys.controller_list);
    l = f_index(oa_sys.name{sys},hvac.components.name);
    hvac.outdoor_air.loop(i,1) = hvac.components.loop(l);%find the loop or zone this component is associated with
    hvac.outdoor_air.control(i,1) = oa_list.controller_name(ln);
    hvac.outdoor_air.min_air_schedule(i,1) = {fix_schedule_name(hvac.outdoor_air.min_air_schedule{i},schedules)};
    hvac.outdoor_air.min_frac_schedule(i,1) = {fix_schedule_name(hvac.outdoor_air.min_frac_schedule{i},schedules)};
    hvac.outdoor_air.max_frac_schedule(i,1) = {fix_schedule_name(hvac.outdoor_air.max_frac_schedule{i},schedules)};
end

%Air routing
[hvac.zone_2_plenum,hvac.zone_2_loop,hvac.plenum_2_loop] = connect_zones_loops(zone_names,hvac.plenum,hvac.mixers,hvac.loop,zone_hvac);

%Managers
manager_list = load_object(objects,'AvailabilityManagerAssignmentList',{''},{{'manager_type','manager_name'}},2);
prop1= {'schedule';'fan_schedule';'control_type';'thermostat_tolerance';'runtime_control_type';'runtime'};
prop2 = {'schedule';};
hvac.managers = load_object(objects,'AvailabilityManager',{'NightCycle';'Scheduled';'ScheduledOn';'ScheduledOff';'NightVentilation';'DifferentialThermostat';'HighTemperatureTurnOff';'HighTemperatureTurnOn';'LowTemperatureTurnOff';'LowTemperatureTurnOn';'HybridVentilation';},{prop1,prop2},0);
hvac.managers.loop = zeros(length(hvac.managers.name),1);
for i = 1:1:length(hvac.managers.name)
    for j = 1:1:length(manager_list.name)
        k = f_index(hvac.managers.name{i},manager_list.manager_name{j});
        if ~isempty(k)
            hvac.managers.loop(i) = f_index(manager_list.name{j},hvac.loop.manager);
            break
        end
    end
    %% can there be managers for zone equipment? (unitary_sys)
    switch hvac.managers.type{i}
        case 'NightCycle'
            hvac.managers.fan_schedule(i) = {fix_schedule_name(hvac.managers.fan_schedule{i},schedules)};
        case 'Scheduled'
            hvac.managers.schedule(i) = {fix_schedule_name(hvac.managers.schedule{i},schedules)};
        otherwise
            disp('new type of manager to load')
    end
end
end%Ends function import_hvac

function design = load_design_size(objects,air_loop)
A = f_index('Sizing:Parameters',objects.type);
design.heat_sizing_factor = str2double(objects.name{A});
design.cooling_sizing_factor = str2double(objects.value{A}{1});
design.averaging_window = str2double(objects.value{A}{2});

A = f_index('SizingPeriod:DesignDay',objects.type);
for k = 1:1:length(A)
    value = objects.value{A(k)};
    if strcmpi(value{3},'WinterDesignDay')
        k2 = 1;
    else
        k2 = 2;
    end
    design.design_day_name(k2,1) = objects.name(A(k));%name
    design.month(k2,1) = str2double(value{1});%month
    design.day(k2,1) = str2double(value{2});%day 
    design.drybulb_temp(k2,1) = str2double(value{4});
    design.drybulb_range(k2,1) = str2double(value{5});
    design.wetbulb_temp(k2,1) = str2double(value{9});
    design.pressure(k2,1) = str2double(value{14});
    design.windspeed(k2,1) = str2double(value{15});
    design.winddir(k2,1) = str2double(value{16});
    design.clearness(k2,1) = str2double(value{25});
end

prop1 = {'sizing_load_type';'outdoor_flow';'heating_max_flow';'T_preheat';'w_preheat';'T_precool';'w_precool';'Tsupply_c';'Tsupply_h';'zone_sum_type';...
    'cool_all_outdoor_air';'heat_all_outdoor_air';'w_supply_c';'w_supply_h';'cooling_supply_flow_rate_method';'cooling_supply_flow_rate';'cooling_supply_flow_rate_area';'cool_frac_autosized';'cooling_supply_flow_rate_per_capacity';...
    'heating_supply_flow_rate';'heating_supply_flow_rate_area';'heat_frac_autosized';'heating_supply_flow_rate_per_capacity';'outdoor_air_method';'cool_capacity_method';'cool_capacity';'cool_capacity_area';'frac_autosize_cooling';...
    'heat_capacity_method';'heat_capacity';'heat_capacity_area';'frac_autosize_heating';'control_method';};
loop_size = load_object(objects,'Sizing',{'System'},{prop1},0);
transfer = {'outdoor_flow';'T_preheat';'w_preheat';'T_precool';'w_precool';'Tsupply_c';'Tsupply_h';'w_supply_c';'w_supply_h';};
    
for i = 1:1:length(loop_size.name)
    k = f_index(loop_size.name{i},air_loop.name);%match sizing to loop #
    for j = 1:1:length(transfer)
        design.(transfer{j})(k,1) = loop_size.(transfer{j})(i);
    end
    if strcmpi(loop_size.cool_all_outdoor_air{i},'No')
        design.all_outdoor_flow_cooling(k,1) = false;
    else
        design.all_outdoor_flow_cooling(k,1) = true;
    end
    if strcmpi(loop_size.heat_all_outdoor_air{i},'No')
        design.all_outdoor_flow_heating(k,1) = false;
    else
        design.all_outdoor_flow_heating(k,1) = true;
    end
    design.control_method(k,1) = loop_size.control_method(i);
end
end%Ends function load_design_size

function components = read_components(branches,mixers,splitters)
%putting components in an order so that they are simulated correctly
components.name = {};
components.type = {};
components.branch = {};
components.branch_number = [];
components.loop = [];
components.inlet = {};
components.outlet = {};
components.bypass = [];
for i = 1:1:length(mixers.name)
    m_in = mixers.inlets{i};
    b = f_index(mixers.outlet{i},branches.inlet);
    for j = 1:1:length(m_in)
        components.name(end+1,1) = mixers.name(i);
        components.type(end+1,1) = {'Mixer'};
        components.inlet(end+1,1) = m_in(j);
        components.outlet(end+1,1) = mixers.outlet(i);
        components.bypass(end+1,1) = false;
        components.branch(end+1,1) = branches.name(b);
        components.branch_number(end+1,1) = b;
        components.loop(end+1,1) = branches.loop(b);
    end
end
j = length(components.name);
for i = 1:1:length(branches.component_type)
    if branches.loop(i)>0
        n = length(branches.component_type{i});
        components.name(j+1:j+n,1) = branches.component_name{i};
        components.type(j+1:j+n,1) = branches.component_type{i};
        components.branch(j+1:j+n,1) = branches.name(i);
        components.branch_number(j+1:j+n,1) = i;
        components.loop(j+1:j+n,1) = branches.loop(i);
        components.inlet(j+1:j+n,1) = branches.component_inlet{i};
        components.outlet(j+1:j+n,1) = branches.component_outlet{i};
        if n ==1 && strcmpi(components.type(j+1),'Pipe:Adiabatic') && any(strcmp(components.outlet(j+1),components.inlet(1:j)) & strcmp('Mixer',components.type(1:j)))
            components.bypass(j+1) = true;
        else
            components.bypass(j+1:j+n,1) = false;
        end
        j = j+n;
    end
end
for i = 1:1:length(splitters.name)
    s_out = splitters.outlets{i};
    b = f_index(splitters.inlet{i},branches.outlet);
    for j = 1:1:length(s_out)
        components.name(end+1,1) = splitters.name(i);
        components.type(end+1,1) = {'Splitter'};
        components.inlet(end+1,1) = splitters.inlet(i);
        components.outlet(end+1,1) = s_out(j);
        components.bypass(end+1,1) = false;
        components.branch(end+1,1) = branches.name(b);
        components.branch_number(end+1,1) = b;
        components.loop(end+1,1) = branches.loop(b);
    end
end
end%ends function read_components

function terminals = load_terminals(objects,splitters,air_loop,e_list)
%splitters connect to demand_inlet
prop1 = {'schedule';'inlet';'max_flow';};
prop2 = {'schedule';'outlet_damper';'inlet';'max_flow';'min_method';'min_frac';'min_flow';'flow_frac_schedule';'reheat_coil_type';'reheat_coil_name';'max_flow_water';'min_flow_water';'outlet';'tolerance';'damper_action';'max_reheat_flow_per_area';'max_reheat_frac';};
terminals = load_object(objects,'AirTerminal:SingleDuct',{'Uncontrolled';'VAV:Reheat'},{prop1,prop2},0);

prop1 = {'outlet';'terminal_type';'terminal';};
air_distribution = load_object(objects,'ZoneHVAC',{'AirDistributionUnit'},{prop1},0);

for i = 1:1:length(terminals.name)
    k = f_index(terminals.name{i},e_list.e_name);
    if isempty(k)
        ad = f_index(terminals.name{i},air_distribution.terminal);
        k = f_index(air_distribution.name{ad},e_list.e_name);
    end
    terminals.zone(i,1) = e_list.zone(k);
    terminals.outlet(i,1) = e_list.outlet(k);
    terminals.return(i,1) = e_list.return(k);
    for j = 1:1:length(splitters.name)
        if any(strcmpi(terminals.inlet{i},splitters.outlets{j}))
            terminals.loop(i,1) = f_index(splitters.inlet{j},air_loop.demand_inlet);
            break
        end
    end
end
end%Ends function load_terminal

function [zone_2_plenum,zone_2_loop,plenum_2_loop] = connect_zones_loops(z_names,plenum,mixers,loop,zone_hvac)
n_p = length(plenum.name);
n_z = length(z_names);
n_pi = 0;
for i = 1:1:n_p
    n_pi = n_pi + (length(plenum.inlets{i})+1);%plenum zone feeds itself (explains +1)
end
n_np = n_z - n_pi;
zone_2_plenum = zeros(n_p+n_np,n_z);
for i = 1:1:n_p
    for j = 1:1:length(plenum.inlets{i})
        h = f_index(plenum.inlets{i}{j},zone_hvac.return);
        z = f_index(zone_hvac.name(h),z_names);
        zone_2_plenum(i,z) = 1;
    end
    z = f_index(plenum.zone{i},z_names);
    zone_2_plenum(i,z) = 1;
end
k = 0;
for i = 1:1:n_z
    if ~any(zone_2_plenum(:,i))
        k = k+1;
        zone_2_plenum(n_p+k,i) = 1;%single zone plenum for those connecting directly to mixer
    end
end

n_m = length(mixers.name);
plenum_2_loop = zeros(n_m,n_p+n_np);
for i = 1:1:n_m %mixer inlets are connected to zone return air or plenum outlet. mixer outlet is connected to air_loop.demand_outlet & return_path.outlet
    k = f_index(mixers.outlet{i},loop.demand_outlet);
    for j = 1:1:length(mixers.inlets{i})
        h = f_index(mixers.inlets{i}{j},zone_hvac.return);
        if isempty(h)
            p = f_index(mixers.inlets{i}{j},plenum.outlet);
            plenum_2_loop(k,p) = 1;
        else
            z = f_index(zone_hvac.name(h),z_names);
            jj = nonzeros((1:n_p+n_np)'.*(zone_2_plenum(:,z)>0));
            plenum_2_loop(k,jj) = 1;
        end
    end
end
zone_2_loop = plenum_2_loop*zone_2_plenum;
end%Ends function connect_zones_loops