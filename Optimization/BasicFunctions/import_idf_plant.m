function [plant,hvac] = import_idf_plant(objects,zones)
j = length(objects.name);
plant.loop = load_plant_loops(objects);

%% need to add condensor loop
A = nonzeros((1:j)'.*strcmp(objects.type,'AirLoopHVAC'));
prop = {'controller';'manager';'supply_flow';'branch_list';'connector_list';'supply_inlet';'demand_outlet';'demand_inlet';'supply_outlet';};
hvac.loop = load_object([],objects.name(A),objects.value(A),prop,0);

A = nonzeros((1:j)'.*strcmp(objects.type,'BranchList'));
branch_list = read_branch_list(objects.name(A),objects.value(A));
A = nonzeros((1:j)'.*strcmp(objects.type,'Connector:Mixer'));
plant.mixers = read_mixers(objects.name(A),objects.value(A));
A = nonzeros((1:j)'.*strcmp(objects.type,'Connector:Splitter'));
plant.splitters = read_splitters(objects.name(A),objects.value(A));

A = nonzeros((1:j)'.*strcmp(objects.type,'CoolingTower:SingleSpeed'));
prop = {'water_inlet';'water_outlet';'water_flow';'air_flow';'fan_power';'UA';'free_convect_flow';'free_convect_factor';'free_convect_UA';'free_convect_UA_factor';'performance_method';'heat_reject_capacity';'nominal_capacity';'free_convect_capacity';'free_convect_size_factor';'des_inlet_dry_bulb';'des_inlet_wet_bulb';'des_approach_temp';'des_range_temp';'basin_heater_capacity';'basin_heater_setpoint';'basin_heater_schedule';'evap_loss_mode';'evap_loss_factor';'drift_loss';'blowdown_mode';'blowdown_ratio';'blowdown_makeup_schedule';'supply_water_tank';'outdoor_air_inlet';'capacity_control';'number_of_cells';'cell_control';'cell_min_water_frac';'cell_max_water_frac';'sizing_factor';};
plant.cooling_tower = load_object([],objects.name(A),objects.value(A),prop,0);

A = nonzeros((1:j)'.*strcmp(objects.type,'Branch'));
[plant.loop,hvac.loop,pl,hv] = sep_branches(branch_list,plant.loop,hvac.loop,objects.name(A),A);

if ~isempty(plant.loop)
    [plant.loop,plant.branches,plant.components] = read_branches(objects.name(pl),objects.value(pl),plant.loop);
    b_o = 0;%branch order #
    branch_order = zeros(length(plant.branches.name),1);
    for p = 1:1:length(plant.loop.name)
        [plant.demand_connect{p},branch_order] = branch_connect(plant.branches,plant.mixers,plant.splitters,plant.loop.demand_outlet{p},plant.loop.demand_branches{p},branch_order,b_o);%% branch connections, working backwards from outlet node)
        b_o = b_o + length(plant.loop.demand_branches{p});
        [plant.supply_connect{p},branch_order] = branch_connect(plant.branches,plant.mixers,plant.splitters,plant.loop.plant_outlet{p},plant.loop.supply_branches{p},branch_order,b_o);%% branch connections, working backwards from outlet node)
        b_o = b_o + length(plant.loop.supply_branches{p});
    end
    [hvac.chiller,hvac.boiler,plant.loop] = load_plant_loop_comp(objects,plant.loop,plant.components);
end

%%hvac portion
[hvac.loop,hvac.branches,hvac.components] = read_branches(objects.name(hv),objects.value(hv),hvac.loop);
hvac.fans = import_fans(objects);
hvac.curves = import_curves(objects);
hvac.design = load_design_size(objects,hvac.loop);
hvac.coils = import_coils(objects,hvac.design,plant);
hvac.outdoor_air = read_outdoor_air(objects,hvac.loop);

prop = {'schedule';'inlet';'outlet';'sensor';'coil_type';'coil_name';};
A = nonzeros((1:j)'.*strcmp(objects.type,'CoilSystem:Cooling:DX'));
hvac.coil_system = load_object([],objects.name(A),objects.value(A),prop,0);
A = nonzeros((1:j)'.*strcmp(objects.type,'AirLoopHVAC:ReturnPath'));
hvac.return_path = read_return_path(objects.name(A),objects.value(A));
A = nonzeros((1:j)'.*strcmp(objects.type,'AirLoopHVAC:ZoneSplitter'));
hvac.splitters = read_splitters(objects.name(A),objects.value(A));
A = nonzeros((1:j)'.*strcmp(objects.type,'AirLoopHVAC:ZoneMixer'));
hvac.mixers = read_mixers(objects.name(A),objects.value(A));
A = nonzeros((1:j)'.*strcmp(objects.type,'AirLoopHVAC:ReturnPlenum'));
hvac.plenum = air_plenum(objects.name(A),objects.value(A));
A = nonzeros((1:j)'.*strcmp(objects.type,'AirLoopHVAC:SupplyPath'));
hvac.supply_path = supply_path(objects.name(A),objects.value(A));

hvac.terminals = load_terminals(objects,zones,hvac.splitters,hvac.loop);

prop = {'schedule';'flow_capacity';'max_power';'fan_power';'standby_power';};
A = nonzeros((1:j)'.*strcmp(objects.type,'Humidifier:Steam:Electric'));
humidifier.type(1:length(A)) = {'Electric'};
hvac.humidifiers = load_object(humidifier,objects.name(A),objects.value(A),prop,0);

hvac.unitary_sys = import_unitary(objects,zones,hvac.loop);

n_p = length(hvac.plenum.name);
n_z = length(zones.name);
n_pi = 0;
for i = 1:1:n_p
    n_pi = n_pi + (length(hvac.plenum.inlets{i})+1);%plenum zone feeds itself (explains +1)
end
n_np = n_z - n_pi;
hvac.zone_2_plenum = zeros(n_p+n_np,n_z);
for i = 1:1:n_p
    for j = 1:1:length(hvac.plenum.inlets{i})
        z = nonzeros((1:length(zones.name))'.*strcmpi(hvac.plenum.inlets{i}{j},zones.return(:,1)));
        hvac.zone_2_plenum(i,z) = 1;
    end
    z = nonzeros((1:length(zones.name))'.*strcmpi(hvac.plenum.zone{i},zones.name));
    hvac.zone_2_plenum(i,z) = 1;
end
k = 0;
for i = 1:1:n_z
    if ~any(hvac.zone_2_plenum(:,i))
        k = k+1;
        hvac.zone_2_plenum(n_p+k,i) = 1;%single zone plenum for those connecting directly to mixer
    end
end

n_m = length(hvac.mixers.name);
hvac.plenum_2_loop = zeros(n_m,n_p+n_np);
for i = 1:1:n_m %mixer inlets are connected to zone return air or plenum outlet. mixer outlet is connected to air_loop.demand_outlet & return_path.outlet
    k = nonzeros((1:length(hvac.loop.name))'.*strcmpi(hvac.mixers.outlet{i},hvac.loop.demand_outlet));
    for j = 1:1:length(hvac.mixers.inlets{i})
        z = nonzeros((1:length(zones.name))'.*strcmpi(hvac.mixers.inlets{i}{j},zones.return(:,1)));
        if isempty(z)
            p = nonzeros((1:length(hvac.plenum.name))'.*strcmpi(hvac.mixers.inlets{i}{j},hvac.plenum.outlet));
            hvac.plenum_2_loop(k,p) = 1;
        else
            jj = nonzeros((1:n_p+n_np)'.*(hvac.zone_2_plenum(:,z)>0));
            hvac.plenum_2_loop(k,jj) = 1;
        end
    end
end
hvac.zone_2_loop = hvac.plenum_2_loop*hvac.zone_2_plenum;
end%Ends function import_idf_plant

function loop = load_plant_loops(objects)
j = length(objects.name);
A = nonzeros((1:j)'.*strcmp(objects.type,'PlantLoop'));
prop = {'fluid';'user_fluid';'operation_scheme';'setpoint_node';'max_temperature';'min_temperature';'max_flow_rate';'min_flow_rate';'loop_volume';'plant_inlet';'plant_outlet';'plant_branch_list';'plant_connector_list';'demand_inlet';'demand_outlet';'demand_branch_list';'demand_connector_list';'load_scheme';};
loop = load_object([],objects.name(A),objects.value(A),prop,0);

B = nonzeros((1:j)'.*strcmp(objects.type,'CondenserLoop'));
loop = load_object(loop,objects.name(B),objects.value(B),prop,length(A));

C = nonzeros((1:j)'.*strcmp(objects.type,'Sizing:Plant'));
for k = 1:1:length(C)
    value = objects.value{C(k)};
    i = nonzeros((1:length(loop.name))'.*strcmpi(objects.name{C(k)},loop.name));
    loop.type(i,1) = value(1);
    loop.exit_temperature(i,1) = str2double(value{2});
    loop.temperature_difference(i,1) = str2double(value{3});
end
end%Ends function load_plant_loops

function [chiller,boiler,loop] = load_plant_loop_comp(objects,loop,components)
j = length(objects.name);
%chillers, boilers
prop = {'capacity';'COP_ref';'chilled_water_temperature_ref';'condenser_water_temperature_ref';'chilled_water_flow_ref';'condenser_water_flow_ref';'capacity_temperature_curve';'electric_in_out_temperature_curve';'electric_in_out_part_load_curve';'min_part_load';'max_part_load';'opt_part_load';'min_unload_ratio';'chilled_water_inlet';'chilled_water_outlet';'condenser_inlet';'condenser_outlet';'frac_power_2_condenser';'min_chilled_water_temperature';'flow_mode_type';'heat_recovery_flow';'heat_recovery_inlet';'heat_recovery_outlet';'sizing_factor';};
A = nonzeros((1:j)'.*strcmp(objects.type,'Chiller:Electric:ReformulatedEIR'));
chiller.type(1:length(A)) = {'ReformulatedEIR'};
chiller = load_object(chiller,objects.name(A),objects.value(A),prop,0);

prop = {'capacity';'COP_ref';'chilled_water_temperature_ref';'condenser_water_temperature_ref';'chilled_water_flow_ref';'condenser_water_flow_ref';'capacity_temperature_curve';'electric_in_out_temperature_curve';'electric_in_out_part_load_curve';'min_part_load';'max_part_load';'opt_part_load';'min_unload_ratio';'chilled_water_inlet';'chilled_water_outlet';'condenser_inlet';'condenser_outlet';'condenser_type';'fan_power_ratio';'frac_power_2_condenser';'min_chilled_water_temperature';'flow_mode_type';'heat_recovery_flow';'heat_recovery_inlet';'heat_recovery_outlet';'sizing_factor';};
B = nonzeros((1:j)'.*strcmp(objects.type,'Chiller:Electric:EIR'));
chiller.type(length(A)+1:length(A)+length(B)) = {'EIR'};
chiller = load_object(chiller,objects.name(B),objects.value(B),prop,length(A));

prop = {'fuel';'capacity';'efficiency';'eff_temperature_curve';'normalizd_eff_curve';'design_water_temperature';'design_flow_rate';'min_part_load';'max_part_load';'optimal_part_load';'water_inlet';'water_outlet';'max_water_temperature';'flow_mode';'parasitic_electric_load';'sizing_factor';};
A = nonzeros((1:j)'.*strcmp(objects.type,'Boiler:HotWater'));
boiler.type(1:length(A)) = {'HotWater'};
boiler = load_object(boiler,objects.name(A),objects.value(A),prop,0);
%update plant loop supply temperatures with chiller.design_water_temperature;
for i = 1:1:length(chiller.type)
    c_num = nonzeros((1:length(components.name))'.*strcmpi(chiller.name{i},components.name));
    l_num = components.loop(c_num);
    loop.cold_water_supply_temperature(l_num(1)) = chiller.chilled_water_temperature_ref(i);
end

for i = 1:1:length(boiler.type)
    c_num = nonzeros((1:length(components.name))'.*strcmpi(boiler.name{i},components.name));
    l_num = components.loop(c_num);
    loop.hot_water_supply_temperature(l_num) = boiler.design_water_temperature(i);
end
end%Ends function load_plant_loop_comp

function [plant,hvac,pl,hv] = sep_branches(branch_list,plant,hvac,names,index)
%%seperate out hvac branches
all_plant_br = [];
if ~isempty(plant)
    for i = 1:1:length(plant.name)
        bl = nonzeros((1:length(branch_list.name))'.*strcmpi(plant.demand_branch_list{i},branch_list.name));
        plant.demand_branches(i,1) = branch_list.branches(bl);
        all_plant_br = [all_plant_br;branch_list.branches{bl}];
        bl = nonzeros((1:length(branch_list.name))'.*strcmpi(plant.plant_branch_list{i},branch_list.name));
        plant.supply_branches(i,1) = branch_list.branches(bl);
        all_plant_br = [all_plant_br;branch_list.branches{bl}];
        plant.branches(i,1) = {[plant.demand_branches{i};plant.supply_branches{i}]};
    end
end
all_hvac_br = [];
if ~isempty(hvac)
    for i = 1:1:length(hvac.name)
        bl = nonzeros((1:length(branch_list.name))'.*strcmpi(hvac.branch_list{i},branch_list.name));
        hvac.branches(i,1) = branch_list.branches(bl);
        all_hvac_br = [all_hvac_br;branch_list.branches{bl}];
    end
end
pl = [];
hv = [];
for i = 1:1:length(names)
    if any(strcmp(names{i},all_hvac_br))
        hv(end+1) = index(i);
    elseif any(strcmp(names{i},all_plant_br))
        pl(end+1) = index(i);
    end
end
end%Ends function sep_branches

function [connect,branch_order] = branch_connect(branches,mixers,splitters,outlet1,loop_branches,branch_order,b_o)%% branch connections, working backwards from outlet node)
n_b = length(loop_branches);
b_num = zeros(n_b,1);
for k = 1:1:n_b
    b_num(k,1) = nonzeros((1:length(branches.name))'.*strcmpi(loop_branches{k},branches.name));
end
    
connect = zeros(n_b,n_b);%matrix row i = from branch, column j = to branch
br_out = nonzeros((1:length(branches.outlet))'.*strcmpi(outlet1,branches.outlet));%all flow assumed to exit from this branch
branch_order(br_out,1) = b_o + n_b;%last demand branch to be simulated
br_in = br_out;
k = 1;
while k<n_b
    br_out = 0*br_in;
    for i = 1:1:length(br_in)
        br = nonzeros((1:length(branches.outlet))'.*strcmpi(branches.inlet{br_in(i)},branches.outlet));%directly connected to another branch
        if ~isempty(br)
            br_out(i) = br;
        else%maybe from splitter
            for j = 1:1:length(splitters.name)
                if any(strcmpi(branches.name{br_in(i)},splitters.outlets{j}))
                    br = nonzeros((1:length(branches.outlet))'.*strcmpi(splitters.inlet{j},branches.name));%directly connected to another branch
                    br_out(i) = br;
                end
            end

        end
    end
    new_br_in = [];
    for i = 1:1:length(br_in)
        if br_out(i)~=0
            if branch_order(br_out(i),1) == 0
                branch_order(br_out(i),1) = b_o + n_b - k;
                k = k+1;
            end
            ii = nonzeros((1:n_b)'.*(b_num == br_out(i)));
            jj = nonzeros((1:n_b)'.*(b_num == br_in(i)));
            connect(ii,jj) = 1;
            new_br_in(end+1) = br_out(i);
        else%branch comes from mixer
            mix = nonzeros((1:length(mixers.outlet))'.*strcmpi(branches.name{br_in(i)},mixers.outlet));
            mix_in = mixers.inlets{mix,1};
            for j = 1:1:length(mix_in)
                br = nonzeros((1:length(branches.name))'.*strcmpi(mix_in{j},branches.name));
                if branch_order(br,1) == 0
                    branch_order(br,1) = b_o + n_b - k;
                    k = k+1;
                end
                ii = nonzeros((1:n_b)'.*(b_num == br(i)));
                jj = nonzeros((1:n_b)'.*(b_num == br_in(i)));
                connect(ii,jj) = 1;
                new_br_in(end+1) = br;
            end
        end
    end
    br_in = new_br_in;
end
end%Ends function branch connect

function unit = import_unitary(objects,zones,loop)
j = length(objects.name);
unit.name = {};
prop = {'schedule';'inlet';'outlet';'fan_type';'fan_name';'max_flow';'heat_coil_type';'heat_coil';'fan_schedule';'fan_without_heating';'max_water';'min_water';'tolerance';'availability';};
A = nonzeros((1:j)'.*strcmp(objects.type,'ZoneHVAC:UnitHeater'));
unit.type(1:length(A)) = {'heater'};
unit = load_object(unit,objects.name(A),objects.value(A),prop,0);

prop = {'schedule';'inlet';'outlet';'fan_schedule';'max_temperature';'cooling_flow';'heating_flow';'no_load_flow';'thermostat';'fan_type';'fan_name';'fan_location';'heat_coil_type';'heat_coil';'cool_coil_type';'cool_coil';'dehumid_control';'reheat_type';'reheat_name';};
B = nonzeros((1:j)'.*strcmp(objects.type,'AirLoopHVAC:UnitaryHeatCool'));
unit.type(length(A)+1:length(A)+length(B)) = {'heat_cool'};
unit = load_object(unit,objects.name(B),objects.value(B),prop,length(A));

n = length(unit.name);
unit.zone = zeros(n,1);
unit.loop = zeros(n,1);
for i = 1:1:n
    k = nonzeros((1:length(loop.name))'.*strcmp(unit.outlet{i},loop.supply_outlet));
    if ~isempty(k)
        unit.loop(i) = k;
    else
        j = 0;
        while isempty(k)
            j = j+1;
            k = nonzeros((1:length(zones.name))'.*strcmp(unit.outlet{i},zones.inlet(:,j)));
        end
        unit.zone(i) = k;
    end
end
end%Ends function import_unitary

function fan = import_fans(objects)
j = length(objects.name);
A = nonzeros((1:j)'.*strcmp(objects.type,'Fan:ConstantVolume'));
prop = {'schedule';'fan_efficiency';'pressure_rise';'flow_rate';'motor_efficiency';'motor_frac';'inlet';'outlet';'end_use';};
fan.type(1:length(A)) = {'ConstantVolume'};
fan = load_object(fan,objects.name(A),objects.value(A),prop,0);

B = nonzeros((1:j)'.*strcmp(objects.type,'Fan:ZoneExhaust'));
prop = {'schedule';'fan_efficiency';'pressure_rise';'flow_rate';'inlet';'outlet';'end_use';};
fan.type(length(A)+1:length(A)+length(B)) = {'ZoneExhaust'};
fan = load_object(fan,objects.name(B),objects.value(B),prop,length(A));

C = nonzeros((1:j)'.*strcmp(objects.type,'Fan:VariableVolume'));
prop = {'schedule';'fan_efficiency';'pressure_rise';'flow_rate';'method';'min_flow_frac';'min_flow_rate';'motor_efficiency';'motor_frac';'a1';'a2';'a3';'a4';'a5';'inlet';'outlet';'end_use';};
fan.type(length(A)+length(B)+1:length(A)+length(B)+length(C)) = {'VariableVolume'};
fan = load_object(fan,objects.name(C),objects.value(C),prop,length(A)+length(B));

D = nonzeros((1:j)'.*strcmp(objects.type,'Fan:OnOff'));
prop = {'schedule';'fan_efficiency';'pressure_rise';'flow_rate';'motor_efficiency';'motor_frac';'inlet';'outlet';'power_speed_curve';'efficiency_speed_curve';'end_use';};
fan.type(length(A)+length(B)+length(C)+1:length(A)+length(B)+length(C)+length(D)) = {'OnOff'};
fan = load_object(fan,objects.name(D),objects.value(D),prop,length(A)+length(B)+length(C));

%convert a few things from cell arrays to numeric arrays
if ~isempty([A;B;C;D;]) && ~isnumeric(fan.flow_rate)
    ffr = fan.flow_rate;
    fan.flow_rate = zeros(length(ffr),1);
    for i = 1:1:length(ffr)
        fan.flow_rate(i) = str2double(ffr{i});
    end
end
end%Ends function import_fans

function coils = import_coils(objects,design,plant)
j = length(objects.name);
A = nonzeros((1:j)'.*strcmp(objects.type,'Coil:Cooling:DX:TwoSpeed'));
if ~isempty(A)
    prop = {'schedule';'rated_capacity';'rated_sensible_heat_ratio';'rated_COP';'rated_air_flow';'static_air_pressure';'inlet';'outlet';'capacity_v_temperature_curve';'capacity_v_flow_curve';'energy_input_v_temperature_curve';'energy_input_v_flow_curve';'part_load_curve';'rated_capacity2';'rated_sensible_heat_ratio2';'rated_COP2';'rated_air_flow2';'capacity_v_temperature_curve2';'energy_input_v_temperature_curve2';'condenser_inlet';};
    coils.cooling.type(1:length(A),1) = {'DX:TwoSpeed'};
    coils.cooling = load_object(coils.cooling,objects.name(A),objects.value(A),prop,0);
    cc_rc2 = coils.cooling.rated_capacity2;
    cc_rf2 = coils.cooling.rated_air_flow2;
    n_cc2 = length(cc_rc2);
    coils.cooling.rated_capacity2 = zeros(n_cc2,1);
    coils.cooling.rated_air_flow2 = zeros(n_cc2,1);
    for i = 1:1:n_cc2
        coils.cooling.rated_capacity2(i) = str2double(cc_rc2{i});
        coils.cooling.rated_air_flow2(i) = str2double(cc_rf2{i});
    end
end

B = nonzeros((1:j)'.*strcmp(objects.type,'Coil:Cooling:DX:SingleSpeed'));
if ~isempty(B)
    prop = {'schedule';'rated_capacity';'rated_sensible_heat_ratio';'rated_COP';'rated_air_flow';'rated_fan_power_per_volume';'inlet';'outlet';'capacity_v_temperature_curve';'capacity_v_flow_curve';'energy_input_v_temperature_curve';'energy_input_v_flow_curve';'part_load_curve';'min_ambient_temp';};
    coils.cooling.type(length(A)+1:length(A)+length(B),1) = {'DX:SingleSpeed'};
    coils.cooling = load_object(coils.cooling,objects.name(B),objects.value(B),prop,length(A));
end
if ~isempty(A) || ~isempty(B)
    hc_rc = coils.cooling.rated_capacity;
    cc_rf = coils.cooling.rated_air_flow;
    cc_shr = coils.cooling.rated_sensible_heat_ratio;

    n_hc = length(hc_rc);
    coils.cooling.rated_capacity = zeros(n_hc,1);
    coils.cooling.rated_air_flow = zeros(n_hc,1);
    coils.cooling.rated_sensible_heat_ratio = zeros(n_hc,1);
    for i = 1:1:n_hc
        coils.cooling.rated_capacity(i) = str2double(hc_rc{i});
        coils.cooling.rated_air_flow(i) = str2double(cc_rf{i});
        coils.cooling.rated_sensible_heat_ratio(i) = str2double(cc_shr{i});
    end
end

C = nonzeros((1:j)'.*strcmp(objects.type,'Coil:Cooling:Water'));
prop = {'schedule';'water_flow';'air_flow';'water_inlet_temperature';'inlet_temperature';'outlet_temperature';'air_inlet_humididty';'air_outlet_humidity';'water_inlet_node';'water_outlet_node';'air_inlet_node';'air_outlet_node';'analysis_type';'heat_exchanger_configuration';};
coils.cooling.type(length(A)+length(B)+1:length(A)+length(B)+length(C),1) = {'Water'};
coils.cooling = load_object(coils.cooling,objects.name(C),objects.value(C),prop,length(A)+length(B));

for i = length(A)+length(B)+1:1:length(A)+length(B)+length(C)% fix some auto-size parameters
    c_num = nonzeros((1:length(plant.components.name))'.*strcmpi(coils.cooling.name{i},plant.components.name));
    if iscell(coils.cooling.water_inlet_temperature(i))
        cw_t = str2double(coils.cooling.water_inlet_temperature{i});
    else
        cw_t = coils.cooling.water_inlet_temperature(i);
    end
    if isnan(cw_t)
        coils.cooling.inlet_water_temperature(i,1) = plant.loop.cold_water_supply_temperature(plant.components.loop(c_num));
    else
        coils.cooling.inlet_water_temperature(i,1) = cw_t;
    end
    if iscell(coils.cooling.outlet_temperature(i))
        ca_t = str2double(coils.cooling.outlet_temperature{i});
    else
        ca_t = coils.cooling.outlet_temperature(i);
    end
    if isnan(ca_t)
        coils.cooling.air_outlet_temperature(i,1) = design.Tsupply_c(plant.components.loop(c_num));
    else
        coils.cooling.air_outlet_temperature(i,1) = ca_t;
    end
    if iscell(coils.cooling.inlet_temperature(i))
        ci_t = str2double(coils.cooling.inlet_temperature{i});
    else
        ci_t = coils.hcooling.inlet_temperature(i);
    end
    coils.cooling.air_inlet_temperature(i,1) = ci_t;
    pl = plant.components.loop(nonzeros((1:length(plant.components.name))'.*strcmpi(coils.cooling.name{i},plant.components.name)));%plant loop
    coils.cooling.outlet_water_temperature(i,1) = coils.cooling.inlet_water_temperature(i,1) + plant.loop.temperature_difference(pl);
    coils.cooling.rated_capacity(i) = nan;
end

D = nonzeros((1:j)'.*strcmp(objects.type,'Coil:Cooling:Water:DetailedGeometry'));
prop = {'schedule';'max_flow';'area_outside';'area_inside';'fin_area';'airflow_area';'coil_depth';'fin_diameter';'tube_id';'tube_od';'tube_conductivity';'fin_conductivity';'fin_spacing';'tube_spacing';'tube_rows';'tubes_per_row';'inlet_water';'outlet_water';'inlet';'outlet';};
coils.cooling.type(length(A)+length(B)+length(C)+1:length(A)+length(B)+length(C)+length(D),1) = {'Water:DetailedGeometry'};
coils.cooling = load_object(coils.cooling,objects.name(D),objects.value(D),prop,length(A)+length(B)+length(C));


for i = length(A)+length(B)+length(C)+1:1:length(A)+length(B)+length(C)+length(D)% fix some auto-size parameters
    c_num = nonzeros((1:length(plant.components.name))'.*strcmpi(coils.cooling.name{i},plant.components.name));
    coils.cooling.inlet_water_temperature(i,1) = plant.loop.cold_water_supply_temperature(plant.components.loop(c_num));
    coils.cooling.air_outlet_temperature(i,1) = design.Tsupply_c(plant.components.loop(c_num));
    pl = plant.components.loop(nonzeros((1:length(plant.components.name))'.*strcmpi(coils.cooling.name{i},plant.components.name)));%plant loop
    coils.cooling.outlet_water_temperature(i,1) = coils.cooling.inlet_water_temperature(i,1) + plant.loop.temperature_difference(pl);
    coils.cooling.rated_capacity(i) = nan;
end

A = nonzeros((1:j)'.*strcmp(objects.type,'Coil:Heating:Fuel'));
prop = {'schedule';'fuel';'efficiency';'capacity';'inlet';'outlet';'temperature_setpoint';};
coils.heating.type(1:length(A),1) = {'Fuel'};
coils.heating = load_object(coils.heating,objects.name(A),objects.value(A),prop,0);

B = nonzeros((1:j)'.*strcmp(objects.type,'Coil:Heating:Water'));
prop = {'schedule';'u_factor_times_area';'max_flow';'inlet_water';'outlet_water';'inlet';'outlet';'input_method';'capacity';'water_inlet_temperature';'inlet_temperature';'water_outlet_temperature';'outlet_temperature';'air_water_convection_ratio';};
coils.heating.type(length(A)+1:length(A)+length(B),1) = {'Water'};
coils.heating = load_object(coils.heating,objects.name(B),objects.value(B),prop,length(A));

C = nonzeros((1:j)'.*strcmp(objects.type,'Coil:Heating:Electric'));
prop = {'schedule';'efficiency';'capacity';'inlet';'outlet';};
coils.heating.type(length(A)+length(B)+1:length(A)+length(B)+length(C),1) = {'Electric'};
coils.heating = load_object(coils.heating,objects.name(C),objects.value(C),prop,length(A)+length(B));

%% fix some auto-size parameters
for i = length(A)+1:1:length(A)+length(B)
    c_num = nonzeros((1:length(plant.components.name))'.*strcmpi(coils.heating.name{i},plant.components.name));
    if iscell(coils.heating.water_inlet_temperature(i))
        hw_t = str2double(coils.heating.water_inlet_temperature{i});
    else
        hw_t = coils.heating.water_inlet_temperature(i);
    end
    if isnan(hw_t)
        coils.heating.inlet_water_temperature(i,1) = plant.loop.hot_water_supply_temperature(plant.components.loop(c_num));
    else
        coils.heating.inlet_water_temperature(i,1) = hw_t;
    end
    if iscell(coils.heating.outlet_temperature(i))
        ha_t = str2double(coils.heating.outlet_temperature{i});
    else
        ha_t = coils.heating.outlet_temperature(i);
    end
    if isnan(ha_t)
        coils.heating.air_outlet_temperature(i,1) = design.Tsupply_h(plant.components.loop(c_num));
    else
        coils.heating.air_outlet_temperature(i,1) = ha_t;
    end
    if iscell(coils.heating.inlet_temperature(i))
        hi_t = str2double(coils.heating.inlet_temperature{i});
    else
        hi_t = coils.heating.inlet_temperature(i);
    end
    coils.heating.air_inlet_temperature(i,1) = hi_t;
    pl = plant.components.loop(nonzeros((1:length(plant.components.name))'.*strcmpi(coils.heating.name{i},plant.components.name)));%plant loop
    coils.heating.outlet_water_temperature(i,1) = coils.heating.inlet_water_temperature(i,1) - plant.loop.temperature_difference(pl);
end

if isfield(coils.heating,'capacity')
    hc_rc = coils.heating.capacity;
    n_hc = length(hc_rc);
    coils.heating.capacity = zeros(n_hc,1);
    for i = 1:1:n_hc
        coils.heating.capacity(i) = str2double(hc_rc{i});
    end
end
end%Ends function import_coils

function outdoor_air = read_outdoor_air(objects,loop)

prop = {'controller_list';'equipment_list';'manager_list';};
A = nonzeros((1:length(objects.name))'.*strcmp(objects.type,'AirLoopHVAC:OutdoorAirSystem'));
oa_sys = load_object([],objects.name(A),objects.value(A),prop,0);

prop = {'controller_type';'controller_name';};
B = nonzeros((1:length(objects.name))'.*strcmp(objects.type,'AirLoopHVAC:ControllerList'));
oa_list = load_object([],objects.name(B),objects.value(B),prop,0);

% prop = {'mixed_air';'outdoor_air';'relief_air';'return_air';};
% C = nonzeros((1:length(objects.name))'.*strcmp(objects.type,'OutdoorAir:Mixer'));
% oa_mix = load_object([],objects.name(C),objects.value(C),prop,0);

% prop = {'mixed_air';'outdoor_air';'relief_air';'return_air';};
% D = nonzeros((1:length(objects.name))'.*strcmp(objects.type,'AirLoopHVAC:OutdoorAirSystem:EquipmentList'));
% oa_equip_list = load_object([],objects.name(C),objects.value(C),prop,0);

prop = {'relief_outlet';'return_air';'mixed_air';'actuator_node';'min_flow';'max_flow';'economizer_control';'economizer_action';'economizer_max_T';'economizer_max_h';'economizer_max_dp';'enthalpy_limit_curve';'economizer_min_T';'lockout';'min_limit_type';'min_air_schedule';'min_frac_schedule';'max_frac_schedule';'mechanical_controller_name';};
E = nonzeros((1:length(objects.name))'.*strcmp(objects.type,'Controller:OutdoorAir'));
oa_control = load_object([],objects.name(E),objects.value(E),prop,0);

for i = 1:1:length(A)
    outdoor_air.name(i,1) = oa_sys.name(i);
    ln = nonzeros((1:length(oa_list.name))'.*strcmp(oa_sys.controller_list{i},oa_list.name));
    con = nonzeros((1:length(oa_control))'.*strcmp(oa_list.controller_name{ln},oa_control.name));
    outdoor_air.control(i,1) = oa_list.controller_name(ln);
    outdoor_air.schedule(i,1) = oa_control.min_air_schedule(con);
    l = [];
    for j = 1:1:length(loop.name)
        l = nonzeros((1:length(loop.component_name{j}))'.*strcmp(outdoor_air.name{i},loop.component_name{j}));
        if ~isempty(l)
            break
        end
    end
    if ~isempty(l)
        outdoor_air.loop(i) = j;
        outdoor_air.zone(i) = 0;
    else
        disp('need condition for outdoor air controller not on loop')
        outdoor_air.loop(i) = 0;
        outdoor_air.zone(i) = z;
    end
end

end%Ends function read_outdoor_air

function obj = load_object(obj,names,values,prop,j1)
for k = 1:1:length(names)
    value = values{k};
    obj.name(j1+k,1) = names(k);
    for i = 1:1:min(length(prop),length(value))
        if ~isfield(obj,prop{i})
            if ~isnan(str2double(value{i}))
                obj.(prop{i}) = zeros(length(names),1);
            else
                obj.(prop{i}) = cell(length(names),1);
            end
        end
        if ~isnumeric(obj.(prop{i}))
            if ~isempty(strfind(prop{i},'schedule'))
                obj.(prop{i})(j1+k,1) = {fix_name(value{i})};
            else
                obj.(prop{i})(j1+k,1) = value(i);
            end
        else
            obj.(prop{i})(j1+k,1) = str2double(value{i});
        end
    end    
end
end%ends function read_object

function list = read_branch_list(names,values)
for k = 1:1:length(names)
    value = values{k};
    list.name(k,1) = names(k);
    list.branches(k,1) = {value};
end
end%ends function read_branch_list

function [loop,branch,components] = read_branches(names,values,loop)
n_l = length(loop.name);
j = 0;
for k = 1:1:length(names)
    value = values{k};
    branch.name(k,1) = names(k);
    branch.pressure_curve(k,1) = {fix_name(value{1})};
    branch.inlet(k,1) = value(1);
    branch.outlet(k,1) = value(end);
    branch.component_type(k,1) = {value(2:4:end)};
    branch.component_name(k,1) = {value(3:4:end)};
    branch.component_inlet(k,1) = {value(4:4:end)};
    branch.component_outlet(k,1) = {value(5:4:end)};
    n = (length(value)-1)/4;
    components.name(j+1:j+n,1) = value(3:4:end);
    components.branch(j+1:j+n,1) = names(k);
    components.branch_number(j+1:j+n,1) = k;
    j = j+n;
end
loop.component_name = cell(n_l,1);
loop.component_type = cell(n_l,1);
n_c = length(components.name);
components.loop = zeros(n_c,1);
for i = 1:1:n_l
    for k = 1:1:length(loop.branches{i})
        br = nonzeros((1:length(branch.name))'.*strcmpi(loop.branches{i}{k},branch.name));
        if isempty(loop.component_name{i})
            loop.component_name(i) = branch.component_name(br); 
            loop.component_type(i) = branch.component_type(br);
        else
            loop.component_name(i) = {[loop.component_name{i};branch.component_name{br}]}; 
            loop.component_type(i) = {[loop.component_type{i};branch.component_type{br}]};
        end 
    end
end

for i = 1:1:n_c%find the plant loop each component is on. Some appear on two loops (cooling and condensor)
    k = 0;
    same_comp = nonzeros((1:n_c)'.*strcmpi(components.name{i},components.name));
    for j = 1:1:n_l
        if any(strcmpi(components.name{i},loop.component_name{j}))
            k = k+1;
            components.loop(same_comp(k),1) = j;
        end
    end
end

end%ends function read_branches

function mix = read_mixers(names,values)
mix = [];
for k = 1:1:length(names)
    value = values{k};
    mix.name(k,1) = names(k);
    mix.outlet(k,1) = value(1);
    mix.inlets(k,1) = {value(2:end)};
end
end%ends function read_mixers

function split = read_splitters(names,values)
split = [];
for k = 1:1:length(names)
    value = values{k};
    split.name(k,1) = names(k);
    split.inlet(k,1) = value(1);
    split.outlets(k,1) = {value(2:end)};
end
end%ends function read_splitters

function ret_path = read_return_path(names,values)
ret_path = [];
for k = 1:1:length(names)
    value = values{k};
    ret_path.name(k,1) = names(k);
    ret_path.outlet(k,1) = value(1);
    ret_path.component_type(k,1) = {value(2:2:end)};
    ret_path.component(k,1) = {value(3:2:end)};
end
end%ends function read_return_path

function supply =supply_path(names,values)
supply = [];
for k = 1:1:length(names)
    value = values{k};
    supply.name(k,1) = names(k);
    supply.inlet(k,1) = value(1);
    supply.component_type(k,1) = {value(2:2:end)};
    supply.component(k,1) = {value(3:2:end)};
end
end%ends function supply_path

function plen = air_plenum(names,values)
plen.name = {};
for k = 1:1:length(names)
    value = values{k};
    plen.name(k,1) = names(k);
    plen.zone(k,1) = value(1);
    plen.zone_node(k,1) = value(2);
    plen.outlet(k,1) = value(3);
    plen.induced_air(k,1) = value(4);
    plen.inlets(k,1) = {value(5:end)};
end
end%ends function air_plenum

function terminals = load_terminals(objects,zones,splitters,air_loop)
%splitters connect to demand_inlet
j = length(objects.name);
A = nonzeros((1:j)'.*strcmp(objects.type,'AirTerminal:SingleDuct:Uncontrolled'));
prop = {'schedule';'inlet';'max_flow';};
terminals.type(1:length(A),1) = {'Uncontrolled'};
terminals = load_object(terminals,objects.name(A),objects.value(A),prop,0);
for i = 1:1:length(A)
    terminals.zone(i,1) = nonzeros((1:length(zones.inlet))'.*strcmpi(terminals.inlet{i},zones.inlet));
    terminals.min_flow(i,1) = 0;
end

% C = nonzeros((1:j)'.*strcmp(objects.type,'ZoneHVAC:AirDistributionUnit'));
% prop = {'outlet';'terminal_type';'terminal';};
% distrib_unit = read_object([],objects.name(C),objects.value(C),prop,0);

B = nonzeros((1:j)'.*strcmp(objects.type,'AirTerminal:SingleDuct:VAV:Reheat'));
prop = {'schedule';'outlet_damper';'inlet';'max_flow';'min_method';'min_frac';'min_flow';'flow_frac_schedule';'reheat_coil_type';'reheat_coil_name';'max_flow_water';'min_flow_water';'outlet';'tolerance';'damper_action';'max_reheat_flow';'max_reheat_frac';};
terminals.type(length(A)+1:length(A)+length(B),1) = {'Reheat'};
terminals = load_object(terminals,objects.name(B),objects.value(B),prop,length(A));
for i = length(A)+1:1:length(A)+length(B)
    terminals.zone(i,1) = nonzeros((1:length(zones.inlet))'.*strcmpi(terminals.outlet{i},zones.inlet));
end

n = length(terminals.name);
for i = 1:1:n
    for j = 1:1:length(splitters.name)
        if any(strcmpi(terminals.inlet{i},splitters.outlets{j}))
            split = j;
        end
    end
    terminals.loop(i,1) = nonzeros((1:length(air_loop.name))'.*strcmpi(splitters.inlet{split},air_loop.demand_inlet));
end
max_flow = nan(n,1);
min_flow = nan(n,1);
min_frac = nan(n,1);
max_water_flow = nan(n,1);
for i = 1:1:n
    if isnumeric(terminals.max_flow(i))
        max_flow(i) = terminals.min_flow(i);
    elseif ~isnan(str2double(terminals.max_flow{i}))
        max_flow(i) = str2double(terminals.max_flow{i});
    end
    if~isempty(B)
        if isfield(terminals,'min_flow') 
            if isnumeric(terminals.min_flow(i))
                min_flow(i) = terminals.min_flow(i);
            elseif ~isnan(str2double(terminals.max_flow{i}))
                min_flow(i) = str2double(terminals.min_flow{i});
            end
        end
        if isfield(terminals,'min_frac') 
            if isnumeric(terminals.min_frac(i))
                min_frac(i) = terminals.min_frac(i);
            elseif ~isnan(str2double(terminals.min_frac{i}))
                min_frac(i) = str2double(terminals.min_frac{i});
            end
        end
        if isfield(terminals,'max_flow_water') 
            if isnumeric(terminals.max_flow_water(i))
                max_water_flow(i) = terminals.max_flow_water(i);
            elseif ~isnan(str2double(terminals.max_flow_water{i}))
                max_water_flow(i) = str2double(terminals.max_flow_water{i});
            end
        end
    end
end
terminals.max_flow = max_flow;
if~isempty(B)
    terminals.min_flow = min_flow;
    terminals.min_frac = min_frac;
    terminals.max_flow_water = max_water_flow;    
end
end%Ends function load_terminal

function curve = import_curves(objects)
j = length(objects.name);
A = nonzeros((1:j)'.*strcmp(objects.type,'Curve:Quadratic'));
prop = {'a0';'a1';'a2';'min_x';'max_x';};
curve.type(1:length(A)) = {'Quadratic'};
curve = load_object(curve,objects.name(A),objects.value(A),prop,0);

B = nonzeros((1:j)'.*strcmp(objects.type,'Curve:Cubic'));
prop = {'a0';'a1';'a2';'a3';'min_x';'max_x';};
curve.type(length(A)+1:length(A)+length(B)) = {'Cubic'};
curve = load_object(curve,objects.name(B),objects.value(B),prop,length(A));

C = nonzeros((1:j)'.*strcmp(objects.type,'Curve:Biquadratic'));
prop = {'a0';'a1';'a2';'b1';'b2';'ab';'min_x';'max_x';'min_y';'max_y';};
curve.type(length(A)+length(B)+1:length(A)+length(B)+length(C)) = {'Biquadratic'};
curve = load_object(curve,objects.name(C),objects.value(C),prop,length(A)+length(B));

D = nonzeros((1:j)'.*strcmp(objects.type,'Curve:Bicubic'));
prop = {'a0';'a1';'a2';'b1';'b2';'ab';'a3';'b3';'aab';'abb';'min_x';'max_x';'min_y';'max_y';};
curve.type(length(A)+length(B)+length(C)+1:length(A)+length(B)+length(C)+length(D)) = {'Bicubic'};
curve = load_object(curve,objects.name(D),objects.value(D),prop,length(A)+length(B)+length(C));
end%Ends function import_curves

function design = load_design_size(objects,air_loop)
j = length(objects.name);
A = nonzeros((1:j)'.*strcmp(objects.type,'Sizing:Parameters'));
design.heat_sizing_factor = str2double(objects.name{A});
design.cooling_sizing_factor = str2double(objects.value{A}{1});
design.averaging_window = str2double(objects.value{A}{2});

A = nonzeros((1:j)'.*strcmp(objects.type,'SizingPeriod:DesignDay'));
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

A = nonzeros((1:j)'.*strcmp(objects.type,'Sizing:System'));
prop = {'sizing_load_type';'outdoor_flow';'heating_max_flow';'T_preheat';'w_preheat';'T_precool';'w_precool';'Tsupply_c';'Tsupply_h';'zone_sum_type';...
    'cool_all_outdoor_air';'heat_all_outdoor_air';'w_supply_c';'w_supply_h';'cooling_supply_flow_rate_method';'cooling_supply_flow_rate';'cooling_supply_flow_rate_area';'cool_frac_autosized';'cooling_supply_flow_rate_per_capacity';...
    'heating_supply_flow_rate';'heating_supply_flow_rate_area';'heat_frac_autosized';'heating_supply_flow_rate_per_capacity';'outdoor_air_method';'cool_capacity_method';'cool_capacity';'cool_capacity_area';'frac_autosize_cooling';...
    'heat_capacity_method';'heat_capacity';'heat_capacity_area';'frac_autosize_heating';'control_method';};
loop_size = load_object([],objects.name(A),objects.value(A),prop,0);
transfer = {'T_preheat';'w_preheat';'T_precool';'w_precool';'Tsupply_c';'Tsupply_h';'w_supply_c';'w_supply_h';};
    
for i = 1:1:length(loop_size.name)
    k = nonzeros((1:length(air_loop.name))'.*strcmpi(loop_size.name{i},air_loop.name));%match sizing to loop #
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

function name = fix_name(name)
name = rmv_spaces(name);
if ~isempty(name)
    name = strrep(strrep(name,' ','_'),',','');
    name = strrep(strrep(name,'-','_'),'/','_');
    name = strrep(strrep(name,':','_'),'#','number');
    name = strrep(strrep(name,'(',''),')','');
    if any(strcmp(name(1),{'1';'2';'3';'4';'5';'6';'7';'8';'9';'0';}))
        name = strcat('a_',name);
    end
    if length(name)>60
        name = name(1:60);
    end
end
end%Ends function fix_name

function name = rmv_spaces(name)
while ~isempty(name) && strcmp(name(1),' ')
    name = name(2:end);
end
while ~isempty(name) && strcmp(name(end),' ')
    name = name(1:end-1);
end
end%Ends function rmv_spaces