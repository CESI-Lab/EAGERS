function building = import_idf(filename)
[objects,ground_temps] = read_objects(filename);
j = length(objects.name);
i = nonzeros((1:j)'.*strcmp(objects.type,'Building'));
building.name = objects.name{i};

i = nonzeros((1:j)'.*strcmp(objects.type,'RunPeriod'));
switch objects.value{i}{5}
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
D1 = datenum([Y,str2double(objects.value{i}{1}),str2double(objects.value{i}{2}),0,0,0]);
D2 = datenum([Y,str2double(objects.value{i}{3}),str2double(objects.value{i}{4}),0,0,0])+1;
building.sim_date = linspace(D1,D2,(D2-D1)*24+1)';

building.location = {};
i = nonzeros((1:j)'.*strcmp(objects.type,'Site:Location'));
building.location.Latitude = str2double(objects.value{i}{1});
building.location.Longitude = str2double(objects.value{i}{2});
building.location.TimeZone = str2double(objects.value{i}{3});
building.location.Elevation = str2double(objects.value{i}{4});
    
A = nonzeros((1:j)'.*strcmp(objects.type,'Schedule:Compact'));
building.schedule = import_idf_schedules(objects.name(A),objects.value(A));

A = nonzeros((1:j)'.*strcmp(objects.type,'RunPeriodControl:SpecialDays'));
building.holidays.name = objects.name(A);
for k = 1:1:length(A)
    building.holidays.start(k,1) = objects.value{A(k)}(1);
    building.holidays.duration(k,1) = str2double(objects.value{A(k)}{2});
    building.holidays.type(k,1) = objects.value{A(k)}(3);
end
building.ground_temperatures = ground_temps;

building.material = import_idf_material(objects);
A = nonzeros((1:j)'.*strcmp(objects.type,'Construction'));
building.construction = import_idf_construction(objects.name(A),objects.value(A),fieldnames(building.material));
building.exterior = import_idf_exterior(objects);
building.surfaces = import_idf_surfaces(objects,building.construction,building.material);
[building.windows,~] = import_window_door(objects,building.construction,building.material);
for w = 1:1:length(building.windows.name)
    k = nonzeros((1:length(building.surfaces.name))'.*strcmpi(building.windows.surf_name{w},building.surfaces.name));
    building.windows.surface(w,1) = k;
    building.windows.zone(w,1) = building.surfaces.zone(k);
    building.windows.height(w,1) =  mean(building.windows.vertices{w,1}(:,3)) - min(building.surfaces.vertices{k}(:,3));%height from floor
    building.surfaces.area(k) = building.surfaces.area(k) - building.windows.area(w,1); %remove window area (or adjust surface resistance based on U-factor)
end
building.surfaces.i_nodes = 5;%interior wall states per layer minimum is 3, both surfaces and center
building.zones = import_idf_zones(objects);
[building.zones,building.surfaces,building.windows] = zone_geometry(building.zones,building.surfaces,building.windows,building.construction,building.material);
building.zones = import_infiltration_mixing_oa(objects,building.zones,building.schedule);
[building.gains,building.zones] = import_idf_gains(objects,building.zones);
building.ctf = build_ctf(building);
[building.plant,building.HVAC] = import_idf_plant(objects,building.zones);

%% fix schedule names
for i = 1:1:length(building.HVAC.fans.name)
    building.HVAC.fans.schedule(i) = {fix_schedule_name(building.HVAC.fans.schedule{i},building.schedule)};
end
%% other component and equipment categories
A = nonzeros((1:j)'.*strcmp(objects.type,'Refrigeration:Case'));
prop = {'schedule';'zone';'rated_ambient_T';'rated_ambient_RH';'capacity_per_length';'latent_heat_ratio';'runtime_fraction';'length';'temperature';'latent_credit_curve_type';'latent_credit_curve_name';'fan_power_per_length';'operating_fan_power_per_length';'standard_lighting_per_unit_length';'installed_standard_lighting_per_unit_length';'light_schedule';'light_to_case';'anti_sweat_heater_per_length';'minimum_anti_sweat_per_length';'anti_sweat_control';'humidity_at_zero_percent';'height';'anti_sweat_heat_to_case';'defrost_power_per_length';'defrost_type';'defrost_schedule';'defrost_drip_down_schedule';'defrost_curve_type';'defrost_curve_name';'return_air_frac';'restock_schedule';'case_credit_fraction_schedule'};
building.cases = load_object([],objects.name(A),objects.value(A),prop,0);
if ~isempty(building.cases)
    zones = building.cases.zone;
    building.cases.zone = zeros(length(A),1);
    for k = 1:1:length(A)%some cases don't have last condition
        z = nonzeros((1:length(building.zones.name))'.*strcmpi(zones{k},building.zones.name));
        building.cases.zone(k) = z;
    end
end

building.walk_ins = [];%need to add this object type

A = nonzeros((1:j)'.*strcmp(objects.type,'Refrigeration:CompressorRack'));
prop = {'heat_reject_location';'design_COP';'COP_curve';'fan_power';'fan_power_temperature_curve';'condensor_type';'water_inlet';'water_outlet';'water_loop_type';'water_condensor_temperature_schedule';'water_flow_rate';'water_max_flow';'water_max_outlet_T';'water_min_inlet_T';'evap_condenser_schedule';'evap_condenser_effectiveness';'evap_condenser_air_flow';'basin_heater_capacity';'basin_heater_setpoint';'water_pump_power';'water_supply_tank';'condenser_air_inlet';'end_use';'case_name';};
building.racks = load_object([],objects.name(A),objects.value(A),prop,0);

A = nonzeros((1:j)'.*strcmp(objects.type,'WaterUse:Equipment'));
building.water_use = import_idf_water_use(objects.name(A),objects.value(A),building.zones.name);
end%Ends function import_idf

function [objects,ground_temps] = read_objects(filename)
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
while i<n_l-1
    i = i+1;
    com = strfind(lines{i},',');
    if ~isempty(com) && isempty(strfind(lines{i},'!'))
        if isempty(lines{i+1})% 
            %skip
            if ~isempty(strfind(lines{i},'Site:GroundTemperature:BuildingSurface,'))
                com = strfind(lines{i},',');
                com(end+1) = strfind(lines{i},';');
                ground_temps = zeros(12,1);
                for k = 1:1:12
                    ground_temps(k) = str2double(lines{i}(com(k)+1:com(k+1)-1));
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

function [gains,zones] = import_idf_gains(objects,zones)
gains.occupancy = import_occupancy(objects,zones);
zones.max_occupancy = zeros(length(zones.name),1);
for i = 1:1:length(gains.occupancy.zone)
    z_i = gains.occupancy.zone(i);
    zones.max_occupancy(z_i) = zones.max_occupancy(z_i) + gains.occupancy.nominal(i)*zones.floor_area(z_i); 
end
A = nonzeros((1:length(objects.name))'.*strcmp(objects.type,'Lights'));
gains.lighting_internal = import_lighting(objects.name(A),objects.value(A),zones,gains.occupancy);
A = nonzeros((1:length(objects.name))'.*strcmp(objects.type,'ElectricEquipment'));
gains.plug_load = import_equipment(objects.name(A),objects.value(A),zones,gains.occupancy);
A = nonzeros((1:length(objects.name))'.*strcmp(objects.type,'GasEquipment'));
gains.gas_load = import_equipment(objects.name(A),objects.value(A),zones,gains.occupancy);
end%Ends function import_idf_gains

function zone = import_infiltration_mixing_oa(objects,zone,schedule)
A = nonzeros((1:length(objects.name))'.*strcmp(objects.type,'ZoneInfiltration:DesignFlowRate'));
zone.infiltration.schedule = cell(length(zone.name),1);
zone.infiltration.nominal = zeros(length(zone.name),1);
for j = 1:1:length(A)
    value = objects.value{A(j)};
    k = nonzeros((1:length(zone.name))'.*strcmpi({fix_name(value{1})},zone.name));%zone associated with infiltration
    zone.infiltration.schedule(k,1) = {fix_name(value{2})};
    switch value{3}
        case 'Flow/Zone'
            zone.infiltration.nominal(k,1) = str2double(value{4});
        case 'Flow/Area'
            zone.infiltration.nominal(k,1) = zone.floor_area(k)*str2double(value{5});
        case 'Flow/ExteriorArea'
            zone.infiltration.nominal(k,1) = zone.exterior_area(k)*str2double(value{6});
        case 'AirChanges/Hour'
            zone.infiltration.nominal(k,1) = zone.volume(k)/3600*str2double(value{7});
    end
end

A = nonzeros((1:length(objects.name))'.*strcmp(objects.type,'ZoneMixing'));
zone.mixing.schedule = cell(length(zone.name),1);
zone.mixing.nominal = cell(length(zone.name),1);
zone.mixing.source = cell(length(zone.name),1);
for j = 1:1:length(A)
    value = objects.value{A(j)};
    k = nonzeros((1:length(zone.name))'.*strcmpi({fix_name(value{1})},zone.name));%zone associated with infiltration
    sched = fix_schedule_name(fix_name(value{2}),schedule);
    zone.mixing.schedule(k) = {[zone.mixing.schedule{k};{sched}]};
    if strcmpi(value{3},'Flow/Zone')
            zone.mixing.nominal(k) = {[zone.mixing.nominal{k}; str2double(value{4})]};
    elseif strcmpi(value{3},'Flow/Area')
            zone.mixing.nominal(k) = {[zone.mixing.nominal{k}; zone.floor_area(k)*str2double(value{5})]};
    elseif strcmpi(value{3},'Flow/Person')
            disp('need to edit module so that mixing is fcn of current occupancy')
            zone.mixing.nominal(k) = {[zone.mixing.nominal{k}; str2double(value{6})]};
    elseif strcmpi(value{3},'AirChanges/Hour')
            zone.mixing.nominal(k) = {[zone.mixing.nominal{k}; zone.volume/3600*str2double(value{7})]};
    else
    disp('need more categories in zoneMixing')
    end
    zone.mixing.source(k) = {[zone.mixing.source{k}; nonzeros((1:length(zone.name))'.*strcmpi({fix_name(value{8})},zone.name))]};
end

A = nonzeros((1:length(objects.name))'.*strcmp(objects.type,'DesignSpecification:OutdoorAir'));
outdoor_air= import_idf_outdoor_air(objects.name(A),objects.value(A));

B = nonzeros((1:length(objects.name))'.*strcmp(objects.type,'ThermostatSetpoint:DualSetpoint'));
n = length(B);
name = objects.name(B);%name
sched = cell(n,2);
for k = 1:1:n
    sched(k,1) = {fix_name(objects.value{B(k)}{1})};%heating
    sched(k,2) = {fix_name(objects.value{B(k)}{2})};%cooling
end

C = nonzeros((1:length(objects.name))'.*strcmp(objects.type,'ZoneControl:Thermostat'));
n = length(C);
zone.thermostat = cell(length(zone.name),2);
for k = 1:1:n
    j = nonzeros((1:length(zone.name))'.*strcmpi(fix_name(objects.value{C(k)}{1}),zone.name));
    k2 = nonzeros((1:length(name))'.*strcmpi(objects.value{C(k)}{4},name));
    zone.thermostat(j,:) = sched(k2,:);
end

D = nonzeros((1:length(objects.name))'.*strcmp(objects.type,'Sizing:Zone'));
z = length(zone.name);
zone.t_supply_c = zeros(z,1);
zone.t_supply_h = zeros(z,1);
zone.humid_supply_c = zeros(z,1);
zone.humid_supply_h = zeros(z,1);
zone.outdoor_flow.method = cell(z,1);
zone.outdoor_flow.value = zeros(z,1);
n = length(D);
for k = 1:1:n
    value = objects.value{D(k)};
    j = nonzeros((1:length(zone.name))'.*strcmpi(fix_name(objects.name{D(k)}),zone.name));
    zone.t_supply_c(j,1) = str2double(value{2});
    zone.t_supply_h(j,1) = str2double(value{5});
    zone.humid_supply_c(j,1) = str2double(value{7});
    zone.humid_supply_h(j,1) = str2double(value{8});
    oa_spec = fix_name(value{9});
    oa = nonzeros((1:length(outdoor_air.name))'.*strcmpi(oa_spec,outdoor_air.name));
    zone.outdoor_flow.method(j,1) = outdoor_air.method(oa);
    zone.outdoor_flow.value(j,1) = outdoor_air.value(oa);   
end
end%Ends function import_infiltration_mixing_oa


function equip = import_equipment(names,values,zone,occ)
equip = [];
for k = 1:1:length(names)
    value = values{k};
    equip.name(k,1) = {fix_name(names{k})};
    z = nonzeros((1:length(zone.name))'.*strcmpi(fix_name(value{1}),zone.name));
    equip.zone(k,1) = z;
    equip.schedule(k,1) = {fix_name(value{2})};
    if strcmpi(value{3},'EquipmentLevel')
        power = str2double(value{4});
    elseif strcmpi(value{3},'Watts/Area')
        power = str2double(value{5})*zone.floor_area(z);
    elseif strcmpi(value{3},'Watts/Person')
        k2 = nonzeros((1:length(occ.name))'.*strcmpi(equip.zone{k},occ.zone));
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

function light = import_lighting(names,values,zone,occ)
for k = 1:1:length(names)
    value = values{k};
    light.name(k,1) = {fix_name(names{k})};
    z = nonzeros((1:length(zone.name))'.*strcmpi(fix_name(value{1}),zone.name));
    light.zone(k,1) = z;
    light.schedule(k,1) = {fix_name(value{2})};
    if strcmpi(value{3},'LightingLevel')
        lighting = str2double(value{4});
    elseif strcmpi(value{3},'Watts/Area')
        lighting = str2double(value{5})*zone.floor_area(z);
    elseif strcmpi(value{3},'Watts/Person')
        k2 = nonzeros((1:length(occ.name))'.*strcmpi(light.zone{k},occ.zone));
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

function occ = import_occupancy(objects,zone)
A = nonzeros((1:length(objects.name))'.*strcmp(objects.type,'People'));
for k = 1:1:length(A)
    value = objects.value{A(k)};
    occ.name(k,1) = {fix_name(objects.name{A(k)})};
    z = nonzeros((1:length(zone.name))'.*strcmpi(fix_name(value{1}),zone.name));
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

function zone = import_idf_zones(objects)
A = nonzeros((1:length(objects.name))'.*strcmp(objects.type,'Zone'));
for k = 1:1:length(A)
    value = objects.value{A(k)};
    zone.name(k,1) = {fix_name(objects.name{A(k)})};
    zone.origin(k,1) = str2double(value{2});
    zone.origin(k,2) = str2double(value{3});
    zone.origin(k,3) = str2double(value{4});
    zone.multiplier(k,1) = str2double(value{6});
end

A = nonzeros((1:length(objects.name))'.*strcmp(objects.type,'NodeList'));%node_list
n_list.name = objects.name(A);
for j = 1:1:length(A)
    n_names = objects.value{A(j)};
    for k = 1:1:length(n_names)
        n_list.node_name(j,k) = n_names(k);
    end
end

A = nonzeros((1:length(objects.name))'.*strcmp(objects.type,'ZoneHVAC:EquipmentConnections'));
prop = {'hvac_equipment';'inlet';'exhaust';'node_name';'return';};
n = length(zone.name);
nn = length(n_list.node_name(1,:));
zone.hvac_equipment = cell(n,nn);
zone.inlet = cell(n,nn);
zone.exhaust = cell(n,nn);
zone.node_name = cell(n,nn);
zone.return = cell(n,nn);
for j = 1:1:length(A)
    val = objects.value{A(j)};
    k = nonzeros((1:n)'.*strcmpi({fix_name(objects.name{A(j)})},zone.name));%zone associated with hvac connections
    for p = 1:1:length(prop)
        k2 = nonzeros((1:length(n_list.name))'.*strcmpi(val{p},n_list.name));
        if isempty(k2)
            zone.(prop{p})(k,1) = objects.value{A(j)}(p);
        else
            zone.(prop{p})(k,:) = n_list.node_name(k2,:);
        end
    end
end
end%Ends function import_idf_zones

function [zone,surf,window] = zone_geometry(zone,surf,window,construct,material)
n_z = length(zone.name);
n_sur = length(surf.name);
n_w = length(window.name);
ceil_area = zeros(n_z,1); 
for i = 1:1:length(surf.name)
    if any(strcmp(surf.type{i},{'Roof';'Ceiling'}))
        k = nonzeros((1:n_z)'.*strcmpi(surf.zone{i},zone.name));
        ceil_area(k) = ceil_area(k) + surf.area(i);
    end
end
floor_height = zeros(n_z,1);
surf.zone_index = zeros(n_sur,1);
zone.floor_area = zeros(n_z,1);
zone.exterior_area = zeros(n_z,1);
zone.ceil_height = zeros(n_z,1);
zone.windows = cell(n_z,1);
zone.surfaces = cell(n_z,1);
for i = 1:1:n_sur
    k = nonzeros((1:n_z)'.*strcmpi(surf.zone{i},zone.name));
    surf.zone_index(i) = k;
    zone.surfaces(k) = {[zone.surfaces{k};i]};
    if ~strcmp(surf.type{i},'InternalMass')
        surf.vertices{i} = surf.vertices{i} + ones(length(surf.vertices{i}(:,1)),1)*zone.origin(k,:);%change 3-d position of surface based on zone "origin"
    end
    if strcmp(surf.type{i},'Floor')
        zone.floor_area(k) = zone.floor_area(k) + surf.area(i);
        floor_height(k) = mean(surf.vertices{i}(:,3));%z-coordinate
    elseif any(strcmp(surf.type{i},{'Roof';'Ceiling'}))
        zone.ceil_height(k) = zone.ceil_height(k) + mean(surf.vertices{i}(:,3))*surf.area(i)/ceil_area(k);%z-coordinate 
    end
    if length(surf.boundary)>=i && strcmp(surf.boundary{i},'Outdoors')
        zone.exterior_area(k) = zone.exterior_area(k) + surf.area(i);
    end
end
zone.ceil_height = max(0,zone.ceil_height-floor_height);
zone.volume = zone.ceil_height.*zone.floor_area;
for i = 1:1:length(window.name)
    k = nonzeros((1:n_z)'.*strcmpi(window.zone{i},zone.name));
    zone.windows(k) = {[zone.windows{k};i]};
    window.zone_index(i,1) = k;
    vert = length(window.vertices{i}(:,1));
    window.vertices{i} = window.vertices{i} + ones(vert,1)*zone.origin(k,:);
    if strcmp(surf.boundary{window.surface(i)},'Outdoors')
        zone.exterior_area(k) = zone.exterior_area(k) + window.area(i);
    end
    surf.area(window.surface(i)) = surf.area(window.surface(i)) - window.area(i);
end

zone.view_factors = zeros(n_sur+n_w);
surf.area_rad = surf.area;%in case some internal mass objects surface area must be 'reduced' to enforce reciprocity (area of internal mass greater than area of surfaces in zone)
for k = 1:1:n_z
    %compute internal view factors in single zone
    sur_k = nonzeros((1:n_sur)'.*(surf.zone_index==k));
    win_k = nonzeros((1:n_w)'.*(window.zone_index==k));
    win_sur = ones(length(win_k),1);
    for j = 1:1:length(win_k)
        win_sur(j) = nonzeros((1:length(sur_k))'.*(sur_k == window.surface(win_k(j))));%of the (6) surfaces in this zone, which has the window
    end
    vf_zone = calc_view_factors(surf.area(sur_k),window.area(win_k),win_sur);
    %%combine all view factors into 1 matrix [surfaces;windows];
    sur = [sur_k;win_k+n_sur];
    for j = 1:1:length(sur)
        zone.view_factors(sur(j),sur) = vf_zone(j,:);
        zone.view_factors(sur,sur(j)) = vf_zone(:,j);
    end
end

for j = 1:1:n_sur
    if ~strcmp(surf.type{j},'InternalMass')
        mat = construct.(fix_name(surf.construct{j}));
        surf.resistance(j,1) = {resistance(mat,material,surf.area(j),surf.i_nodes)};
        surf.capacitance(j,1) = {capacitance(mat,material,surf.area(j),surf.resistance(j),surf.i_nodes)};
    else
        surf.capacitance(j,1) = {capacitance(construct.(fix_name(surf.construct{j})),material,surf.area(j),[],surf.i_nodes)};
    end
end
end%Ends function zone_geometry

function outdoor_air = import_idf_outdoor_air(names,values)
for k = 1:1:length(names)
    outdoor_air.name(k,1) = {fix_name(names{k})};
    outdoor_air.method(k,1) = values{k}(1);
    if strcmpi(values{k}{1},'Flow/Person')
            j = 2;
    elseif strcmpi(values{k}{1},'Flow/Area')
            j = 3;
    elseif strcmpi(values{k}{1},'Flow/Zone')
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

function obj = load_object(obj,names,values,prop,j1)
for k = 1:1:length(names)
    value = values{k};
    obj.name(j1+k) = names(k);
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
end%ends function load_object

function surf = import_idf_surfaces(objects,construct,material)
j = length(objects.name);
A = nonzeros((1:j)'.*strcmp(objects.type,'BuildingSurface:Detailed'));
A2 = nonzeros((1:j)'.*strcmp(objects.type,'InternalMass'));
prop1 = {'type';'construct';'zone';'boundary';'object';};%'sun';'wind';'view2ground';'vertices'};
for j = 1:1:length(A)
    value = objects.value{A(j)};
    surf.name(j,1) = {fix_name(objects.name{A(j)})};
    for p = 1:1:length(prop1)
        surf.(prop1{p})(j,1) = {fix_name(value{p})};
    end
    surf.vertices(j,1) = {surface_vertices(value(10:9+str2double(value{9})))};
    [surf.area(j,1),surf.normal(j,:)] = find_area(surf.vertices{j});
    mat = construct.(fix_name(surf.construct{j}));
    surf.roughness(j,1) = {material.(mat{1}).roughness{1}};%outside surface material roughness
    [surf.absorptance.solar(j,:),surf.absorptance.thermal(j,:),surf.absorptance.visible(j,:)] = absorptance(mat,material);
end
surf.roughness_factor(strcmpi(surf.roughness,'VeryRough'),1) = 2.17;
surf.roughness_factor(strcmpi(surf.roughness,'Rough'),1) = 1.67;
surf.roughness_factor(strcmpi(surf.roughness,'MediumRough'),1) = 1.52;
surf.roughness_factor(strcmpi(surf.roughness,'MediumSmooth'),1) = 1.13;
surf.roughness_factor(strcmpi(surf.roughness,'Smooth'),1) = 1.11;

for k = 1:1:length(A2)
    value = objects.value{A2(k)};
    surf.name(j+k,1) = {fix_name(objects.name{A2(k)})};
    surf.type(j+k,1) = {'InternalMass'};
    surf.construct(j+k,1) = {fix_name(value{1})};
    surf.zone(j+k,1) = {fix_name(value{2})};
    surf.area(j+k,1) = str2double(value{3});
    surf.normal(j+k,:) = [0 0 1]; %facing upward
    surf.vertices(j+k) = {''};
    surf.boundary(j+k,1) = {''};
    surf.object(j+k,1) = {''};    
end
end%ends function import_idf_surfaces

function [window,door] = import_window_door(objects,construct,material)
%%frames
A = nonzeros((1:length(objects.name))'.*strcmp(objects.type,'WindowProperty:FrameAndDivider'));
frame_names = cell(length(A),1);
frames = cell(length(A),1);
prop = {'width';'projection_out';'projection_in';'conductance';'conductance_ratio_edge_center';'solar_absoprtance';'visible_absorptance';'emissivity';'divider_type';'divider_width';'horizontal_dividers';'vertical_dividers';'divider_projection_out';'divider_projection_in';'divider_conductance';'divider_conductance_ratio';'divider_solar_absoprtnace';'divider_visible absorptance';'divider_emissivity'};
for i = 1:1:length(A)
    frame_names(i) = objects.name(A(i));
    for j = 1:1:length(prop)
        fr.(prop{j}) = str2double(objects.value{A(i)}{j});
    end
    frames(i) = {fr};
end

%%Windows & doors
B = nonzeros((1:length(objects.name))'.*strcmp(objects.type,'FenestrationSurface:Detailed'));
w = 0;
d = 0;
door = [];
for j = 1:1:length(B)
    value = objects.value{B(j)};
    layers = construct.(fix_name(value{2}));
    if length(layers) == 1
        mat = material.(layers{1});
        if strcmp(mat.type,'simple_glazing') || strcmp(mat.type,'glazing')
            w = w+1;
            window.name(w,1) = objects.name(B(j));
            window.type(w,1) = value(1);
            window.construct(w,1) = value(2);
            window.surf_name(w,1) = {fix_name(value{3})};
            window.frame_name(w,1) = value(7);
            window.vertices(w,1) = {surface_vertices(value(10:9+str2double(value{9})))};
            [window.area(w,1),window.normal(w,:)] = find_area(window.vertices{w,1});
            k = nonzeros((1:length(frames))'.*strcmp(window.frame_name{w},frame_names));
            if ~isempty(k)
                window.frame = frames{k};
            end
            if strcmp(mat.type,'glazing')
                %convert the properties below into correct values for model
    %             prop = {'optical_data_type';'spectral_data_set';'thickness';'transmittance_normal';'reflectance_normal_front';'reflectance_normal_back';'transmittance_normal_visible';'reflectance_normal_visible_front';'reflectance_normal_visible_back';'transmittance_normal_infrared';'emissivity_infrared_front';'emissivity_infrared_back';'conductivity';'dirt_correction';'solar_diffusing';};
            elseif strcmp(mat.type,'simple_glazing')
                window = simple_glazing(window,mat,w);
            end
        else
            %ignore doors for now (keeps properties of surface that door is a part of
            d = d+1;
        end
    else
        %need double pane window capability
    end
end
end%Ends function import_window_door

function vertices = surface_vertices(vert)
vertices = zeros(length(vert),3);
for j = 1:1:length(vert)
    com = strfind(vert{j},',');
    vertices(j,1) = str2double(rmv_spaces(vert{j}(1:com(1)-1)));
    vertices(j,2) = str2double(rmv_spaces(vert{j}(com(1)+1:com(2)-1)));
    vertices(j,3) = str2double(rmv_spaces(vert{j}(com(2)+1:end)));
end
end%Ends function surface_vertices

function material = import_idf_material(objects)
j = length(objects.name);
A = nonzeros((1:j)'.*strcmp(objects.type,'Material'));
prop = {'roughness';'thickness';'conductivity';'density';'specheat';'thermal_absorptance';'solar_absorptance';'visible_absorptance';};
material = mat_prop([],objects.name(A),objects.value(A),prop,'standard');

B = nonzeros((1:j)'.*strcmp(objects.type,'Material:NoMass'));
prop = {'roughness';'thermal_resistance';'thermal_absorptance';'solar_absorptance';};%'visible_absorptance';};
material = mat_prop(material,objects.name(B),objects.value(B),prop,'no mass');

C = nonzeros((1:j)'.*strcmp(objects.type,'WindowMaterial:SimpleGlazingSystem'));
prop = {'u_factor';'solar_heat_gain';'visible_transmittance';};
material = mat_prop(material,objects.name(C),objects.value(C),prop,'simple_glazing');

D = nonzeros((1:j)'.*strcmp(objects.type,'WindowMaterial:Glazing'));
prop = {'optical_data_type';'spectral_data_set';'thickness';'transmittance_normal';'reflectance_normal_front';'reflectance_normal_back';'transmittance_normal_visible';'reflectance_normal_visible_front';'reflectance_normal_visible_back';'transmittance_normal_infrared';'emissivity_infrared_front';'emissivity_infrared_back';'conductivity';'dirt_correction';'solar_diffusing';};
material = mat_prop(material,objects.name(D),objects.value(D),prop,'glazing');

E = nonzeros((1:j)'.*strcmp(objects.type,'WindowMaterial:Gas'));
prop = {'gas';'thickness';};
material = mat_prop(material,objects.name(E),objects.value(E),prop,'gas');
end%ends function import_idf_material

function material = mat_prop(material,names,value,prop,type)
for k = 1:1:length(value)
    name = fix_name(names{k});
    material.(name).type = type;
    for j = 1:1:length(prop)
        if isempty(value{k}(j))
            material.(name).(prop{j}) = nan;
        else
            num = str2double(value{k}{j});
            if isnan(num)
                if strcmp(value{k}{j},'N/A')
                    material.(name).(prop{j}) = num;
                else
                    material.(name).(prop{j}) = {fix_name(value{k}{j})};
                end
            else
                material.(name).(prop{j}) = num;
            end
        end
    end
end
end%ends function mat_prop

function construct = import_idf_construction(names,value,materials)
for k = 1:1:length(value)
    name = fix_name(names{k});
    i = 1;
    layers = {};
    while i<=length(value{k})
        lay = fix_name(value{k}{i});
        if ~any(strcmp(lay,materials))
            m = nonzeros((1:length(materials))'.*strcmpi(lay,materials));
            lay = materials{m};
        end
        layers(end+1) = {lay};
        i = i+1;
    end
    construct.(name) = layers;
end
end%ends function import_idf_construction

function water_use = import_idf_water_use(names,values,zone_names)
prop = {'end_use_category';'peak_flow';'flow_schedule';'target_temperature_schedule';'hot_supply_temperature_schedule';'cold_supply_temperature_schedule';'zone';'sensible_frac_schedule';'latent_frac_schedule'};
water_use = [];
for k = 1:1:length(names)
    value = values{k};
    water_use.name(k,1) = names(k);
    for i = 1:1:length(value)
        if ~isfield(water_use,prop{i})
            if ~isnan(str2double(value{i}))
                water_use.(prop{i}) = zeros(length(names),1);
            else
                water_use.(prop{i}) = cell(length(names),1);
            end
        end
        if ~isnumeric(water_use.(prop{i}))
            water_use.(prop{i})(k,1) = value(i);
        else
            water_use.(prop{i})(k,1) = str2double(value{i});
        end
    end
end
if ~isempty(water_use)
    w_zones = water_use.zone;
    water_use.zone = zeros(length(names),1);
    for j = 1:1:length(names)%some cases don't have last condition
        z = nonzeros((1:length(zone_names))'.*strcmpi(fix_name(w_zones{j}),zone_names));
        water_use.zone(j) = z;
        water_use.flow_schedule(j) = {fix_name(water_use.flow_schedule{j})};
        water_use.target_temperature_schedule(j) = {fix_name(water_use.target_temperature_schedule{j})};
        water_use.hot_supply_temperature_schedule(j) = {fix_name(water_use.hot_supply_temperature_schedule{j})};
        water_use.cold_supply_temperature_schedule(j) = {fix_name(water_use.cold_supply_temperature_schedule{j})};
        water_use.sensible_frac_schedule(j) = {fix_name(water_use.sensible_frac_schedule{j})};
        water_use.latent_frac_schedule(j) = {fix_name(water_use.latent_frac_schedule{j})};
    end
end
end%Ends function import_idf_water_use

function ext = import_idf_exterior(objects)
A = nonzeros((1:length(objects.name))'.*strcmp(objects.type,'Exterior:Lights'));
i = 0;
for k = 1:1:length(A)
    i = i+1;
    ext.name(i,1) = {fix_name(objects.name{A(i)})};
    ext.type(i,1) = {'Lights'};
    value = objects.value{A(k)};
    parameter = objects.parameter{A(k)};
    for j = 1:1:length(value)
        if strcmpi(parameter{j},'Fuel Use Type')
            ext.fuel(i,1) = value(j);
        elseif strcmpi(parameter{j},'Schedule Name')
            ext.schedule{i,1} = value(j);
        elseif strcmpi(parameter{j},'Design Level {W}')
            ext.nominal(i,1) = str2double(value{j});
        elseif strcmpi(parameter{j},'End-Use Subcategory')
            ext.category(i,1) = str2double(value{j});
        elseif strcmpi(parameter{j},'Control Option')
            ext.control(i,1) = str2double(value{j});
        else
            disp('need more categories in exterior:Lights')
        end
    end
end
A = nonzeros((1:length(objects.name))'.*strcmp(objects.type,'Exterior:FuelEquipment'));
for k = 1:1:length(A)
    i = i+1;
    ext.name(i,1) = {fix_name(objects.name{A(k)})};
    ext.type(i,1) = {'Lights'};
    value = objects.value{A(k)};
    parameter = objects.parameter{A(k)};
    for j = 1:1:length(value)
        if strcmpi(parameter{j},'Fuel Use Type')
            ext.fuel(i,1) = value(j);
        elseif strcmpi(parameter{j},'Schedule Name')
            ext.schedule{i,1} = value(j);
        elseif strcmpi(parameter{j},'Design Level {W}')
            ext.nominal(i,1) = str2double(value{j});
        elseif strcmpi(parameter{j},'End-Use Subcategory')
            ext.category(i,1) = str2double(value{j});
        else
            disp('need more categories in exterior:fuelEquipment')
        end
    end
end
end%ends function import_idf_exterior

function sched = import_idf_schedules(names,value)
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
            sched.(name).seasons(end+1) = sum(dy(1:(month-1)))*24 + 24*day; 8760;
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
            mat = zeros(0,2);
            while i<=n_v && strcmpi(value{k}{i}(1:5),'Until')
                col = strfind(value{k}{i},':');
                com = strfind(value{k}{i},',');
                mat(end+1,1) = str2double(value{k}{i}(col(1)+1:col(2)-1)) + str2double(value{k}{i}(col(2)+1:com(1)-1))/60;
                mat(end,2) = str2double(value{k}{i}(com(1)+1:end));
                i = i+1;
            end
            for j = 1:1:length(day)
                sched.(name).(day{j}){season} = [0,mat(1,2);mat;];
            end
        else
            i = i+1;%skip miscelaneous line
        end
    end
end
end%ends function import_idf_schedules

function R = resistance(layers,material,A,i_nodes)
n = length(layers);
R = zeros(n*(i_nodes-1),1);
j = 0;
for i = 1:1:n
    mat = material.(layers{i});
    if isfield(mat,'thickness')
        R(j+1:j+i_nodes-1) = mat.thickness/mat.conductivity/A/(i_nodes-1);% K/W 
    elseif isfield(mat,'thermal_resistance')
        R(j+1:j+i_nodes-1) = mat.thermal_resistance/A/(i_nodes-1);% K/W 
    end
    j = j + i_nodes - 1;
end
%average
R = sum(R)/(i_nodes-1)*ones(1,i_nodes-1);%resistances in series add up
end%Ends function resistance

function C = capacitance(layers,material,A,R,i_nodes)
n = length(layers);
C = zeros(n*(i_nodes-1)+1,1);
j = 0;
for i = 1:1:n
    mat = material.(layers{i});
    if isfield(mat,'thickness')
        cap = mat.density*mat.specheat*mat.thickness*A;%capacitance in J/K
    else %air wall, find thickness that matches R and use that for capacitance
        k_air = 2.5e-2; %W/m*K
        rho_air = 1.225;%kg/m^3
        Cp_air = 1006;%Specific Heat {J/kg-K}
        t = sum(R)*A*k_air;%thickness in m
        cap = t*A*rho_air*Cp_air; %capacitance in J/K
    end
    C(j+1) = C(j+1) + cap/(i_nodes-1)/2;
    C(j+2:j+i_nodes-1) = cap/(i_nodes-1);
    C(j+i_nodes) = cap/(i_nodes-1)/2;
    j = j + i_nodes-1;
end
%average
C = sum(C)/i_nodes*ones(1,i_nodes);
end%Ends function capacitance

function [solar,thermal,visible] = absorptance(layers,material)
mat = material.(layers{end});%inside surface material
solar = zeros(1,2);
thermal = zeros(1,2);
visible = zeros(1,2);
db = ones(1,2);
if isfield(mat,'solar_absorptance')
    solar = db*mat.solar_absorptance;
end
if isfield(mat,'thermal_absorptance')
    thermal = db*mat.thermal_absorptance;
end
if isfield(mat,'visible_absorptance')
    visible = db*mat.visible_absorptance;
end
if length(layers)>1%record properties of outside layer if different
    mat = material.(layers{1});%inside surface material
    if isfield(mat,'solar_absorptance')
        solar(1,2) = mat.solar_absorptance;
    end
    if isfield(mat,'thermal_absorptance')
        thermal(1,2) = mat.thermal_absorptance;
    end
    if isfield(mat,'visible_absorptance')
        visible(1,2) = mat.visible_absorptance;
    end
end
end%Ends function absorptance

function [area,normal] = find_area(vertices)
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
normal_v = cross((P2-P1),(P3-P2));
a = normal_v(1);%x-value of normal vector 
b = normal_v(2);%y-value of normal vector 
c = normal_v(3);%z-value of normal vector 
mag = (a^2 + b^2 + c^2)^.5;%magnitude of normal vector
normal = normal_v/mag;
if a~=0 || b~=0 %not in X-Y plane, need to rotate
    if dot([a,b,0],[1,0,0]) == 0 %single rotation about x-axis
        theta_z = sign(normal(2))*acos(dot(normal,[0,0,1]));%angle between normal and z-axis (rotate around y-axis)
        R_x = [1 0 0; 0 cos(theta_z) -sin(theta_z); 0 sin(theta_z) cos(theta_z);];
        for k = 1:1:length(vertices(:,1))
            vertices(k,:) = (R_x*vertices(k,:)')';
        end
    else
        theta_x = -sign(normal(2))*acos(dot([a,b,0],[1,0,0])/(a^2+b^2)^.5);%angle between shadow and x-axis (rotate around z-axis)
        R_z = [cos(theta_x) -sin(theta_x) 0; sin(theta_x) cos(theta_x) 0; 0 0 1];
        normal2 = (R_z*normal')';
        theta_z = -sign(normal(1))*acos(dot(normal2,[0,0,1]));%angle between normal and z-axis (rotate around y-axis)
        R_y = [cos(theta_z) 0 sin(theta_z); 0 1 0; -sin(theta_z) 0 cos(theta_z);];
        for k = 1:1:length(vertices(:,1))
            vertices(k,:) = (R_y*R_z*vertices(k,:)')';
        end
    end
end
area = polyarea(vertices(:,1),vertices(:,2));
end%Ends function find_area

function name = fix_schedule_name(name,schedule)
f = fieldnames(schedule);
if ~any(strcmp(name,f))
    k = strcmpi(name,f);
    name = f{k};
end
end%Ends function fix_schedule_name

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

function [vf,area] = calc_view_factors(area_s,area_w,win_surf)
n = length(area_s);
w = length(area_w);
vf = zeros(n+w,n+w);
area = [area_s;area_w];
net_area = sum(area);
for i = 1:1:n+w
    seen_area = net_area - area(i) - sum(area_w.*(win_surf==i));
    if i>n
        seen_area = seen_area - area_s(win_surf(i-n));
    end
    for j = 1:1:n+w
        if (j>n && win_surf(j-n)==i) || (i>n && win_surf(i-n)==j)
            %do nothing window and surface dont see each other
        elseif i~=j
            vf(i,j) = area(j)/seen_area;
        end
    end
end
%check reciprocity!: 
area_i = area*ones(1,n+w);
perc_error = ones(n+w,n+w);
while any(any(abs(perc_error)>1e-4))
    perc_error = (vf.*area_i - (vf.*area_i)')./(vf.*area_i);
    perc_error(isnan(perc_error)) = 0;
    vf = vf.*(1-perc_error);
    vf = vf.*((1./sum(vf,2))*ones(1,n+w));%ensure all add up to 1
end
end%Ends function calc_view_factors

function window = simple_glazing(window,mat,w)
%conversion to single pane, see section 7.7, page 285 of reference manual
window.solar_heat_gain(w,1) = mat.solar_heat_gain;
if isempty(mat.visible_transmittance{1})
    window.visible_transmittance(w,1) =nan;
else
    window.visible_transmittance(w,1) = str2double(mat.visible_transmittance{1});
end
if mat.u_factor<5.85
    window.interior_glazing_resistance_w(w,1) = 1/(0.359073*log(mat.u_factor) + 6.949915);
else
    window.interior_glazing_resistance_w(w,1) = 1/(1.788041*mat.u_factor - 2.886625);
end
window.exterior_glazing_resistance_w(w,1) = 1/(0.025342*mat.u_factor + 29.163853);
window.bare_glass_resistance(w,1) = 1/mat.u_factor - window.interior_glazing_resistance_w(w,1) - window.exterior_glazing_resistance_w(w,1);
if window.bare_glass_resistance(w,1)>7
    window.thickness(w,1) = 0.002;
else
    window.thickness(w,1) = 0.05914 - 0.00714/window.bare_glass_resistance(w,1);
end
window.thermal_conductivity(w,1) = window.thickness(w,1)/window.bare_glass_resistance(w,1);%%%

if mat.solar_heat_gain<0.7206
    a =  0.939998*mat.solar_heat_gain^2 + 0.20332*mat.solar_heat_gain;
elseif mat.solar_heat_gain>=0.7206
    a = 1.30415*mat.solar_heat_gain - 0.30515;
end
if mat.solar_heat_gain<=0.15
    b = 0.41040*mat.solar_heat_gain;
elseif mat.solar_heat_gain>0.15
    b = 0.085775*mat.solar_heat_gain^2 + 0.963954*mat.solar_heat_gain - 0.084958;
end
if mat.u_factor>4.5
    window.solar_transmittance(w,1) = a;
elseif mat.u_factor<3.4 
    window.solar_transmittance(w,1) = b;
else
    window.solar_transmittance(w,1) = (mat.u_factor - 3.4)/1.1*(a-b)+b;
end
if isnan(window.visible_transmittance(w,1))
    window.visible_transmittance(w,1) = window.solar_transmittance(w,1);
end
SHGC_Tsol = mat.solar_heat_gain - window.solar_transmittance(w,1);
a = 1/(29.436546*(SHGC_Tsol)^3 - 21.943415*(SHGC_Tsol)^2 + 9.945872*(SHGC_Tsol) + 7.426151);
b = 1/(199.8208128*(SHGC_Tsol)^3 - 90.639733*(SHGC_Tsol)^2 + 19.737055*(SHGC_Tsol) + 6.766575);
c = 1/(2.225824*(SHGC_Tsol) + 20.57708);
d = 1/(5.763355*(SHGC_Tsol) + 20.541528);
if mat.u_factor>4.5
    window.interior_glazing_resistance_s(w,1) = a;
    window.exterior_glazing_resistance_s(w,1) = c;
elseif mat.u_factor<3.4 
    window.interior_glazing_resistance_s(w,1) = b;
    window.exterior_glazing_resistance_s(w,1) = d;
else
    window.interior_glazing_resistance_s(w,1) = (mat.u_factor - 3.4)/1.1*(a-b)+b;
    window.exterior_glazing_resistance_s(w,1) = (mat.u_factor - 3.4)/1.1*(c-d)+d;
end
inward_fraction = (window.exterior_glazing_resistance_s(w,1) + 0.5*window.bare_glass_resistance(w,1))/(window.exterior_glazing_resistance_s(w,1) + window.bare_glass_resistance(w,1) + window.interior_glazing_resistance_s(w,1));
window.solar_reflectance(w,1) = 1 - window.solar_transmittance(w,1) - SHGC_Tsol/inward_fraction;
window.thermal_absorptance(w,1) = 0.84;%%%
window.emittance_front(w,1) = 0.84;
window.emittance_back(w,1) = 0.84;
window.long_wave_transmittance(w,1) = 0;
window.visible_reflectance_back(w,1) = -0.7409*window.visible_transmittance(w,1)^3 + 1.6531*window.visible_transmittance(w,1)^2 - 1.2299*window.visible_transmittance(w,1) + 0.4547;
window.visible_reflectance_front(w,1) = -0.0622*window.visible_transmittance(w,1)^3 + 0.4277*window.visible_transmittance(w,1)^2 - 0.4169*window.visible_transmittance(w,1) + 0.2399;
[window.transmittance(w,1:5),window.reflectance(w,1:5)] = transmit_reflect_coef(mat.solar_heat_gain,mat.u_factor);%%%%
end%Ends function simple_glazing

function [transmittance,reflectance] = transmit_reflect_coef(SHGC,U)
%% Is not Used
%Details in engineering reference pg 287 and http://gaia.lbl.gov/btech/papers/2804.pdf
T = [1.470E-02 1.486E+00 -3.852E+00 3.355E+00 -1.474E-03;
     5.546E-01 3.563E-02 -2.416E+00 2.831E+00 -2.037E-03;
     7.709E-01 -6.383E-01 -1.576E+00 2.448E+00 -2.042E-03;
     3.462E-01 3.963E-01 -2.582E+00 2.845E+00 -2.804E-04;
     2.883E+00 -5.873E+00 2.489E+00 1.510E+00 -2.577E-03;
     3.025E+00 -6.366E+00 3.137E+00 1.213E+00 -1.367E-03;
     3.229E+00 -6.844E+00 3.535E+00 1.088E+00 -2.891E-03;
     3.334E+00 -7.131E+00 3.829E+00 9.766E-01 -2.952E-03;
     3.146E+00 -6.855E+00 3.931E+00 7.860E-01 -2.934E-03;
     3.744E+00 -8.836E+00 6.018E+00 8.407E-02 4.825E-04;];

R = [1.632E+01 -5.782E+01 7.924E+01 -5.008E+01 1.334E+01;
     4.048E+01 -1.193E+02 1.348E+02 -7.097E+01 1.611E+01;
     5.749E+01 -1.645E+02 1.780E+02 -8.875E+01 1.884E+01;
     5.714E+00 -1.667E+01 1.863E+01 -9.756E+00 3.074E+00;
     -5.488E-01 -6.498E+00 2.120E+01 -2.097E+01 7.814E+00;
     4.290E+00 -1.267E+01 1.466E+01 -8.153E+00 2.871E+00;
     2.174E+01 -6.444E+01 7.489E+01 -4.179E+01 1.062E+01;
     4.341E+00 -1.280E+01 1.478E+01 -8.203E+00 2.879E+00;
     4.136E+01 -1.178E+02 1.276E+02 -6.437E+01 1.426E+01;
     4.490E+00 -1.266E+01 1.397E+01 -7.501E+00 2.693E+00;]; 

T_fghi = .25*T(6,:) + .25*T(7,:) + .25*T(8,:) + .25*T(9,:); 
R_fghi = .25*R(6,:) + .25*R(7,:) + .25*R(8,:) + .25*R(9,:);
T_bdcd = .25*T(2,:) + .25*T(3,:) + .5*T(4,:);
R_bdcd = .25*R(2,:) + .25*R(3,:) + .5*R(4,:);
T_fh = 0.5*T(6,:) + 0.5*T(8,:);
R_fh = 0.5*R(6,:) + 0.5*R(8,:);
if SHGC>.65 
    if U>4.54
        transmittance = T(1,:);
        reflectance = R(1,:);
    elseif U>3.41
        a = (U-3.41)/1.13;
        transmittance = a*T(1,:) + (1-a)*T(5,:);    
        reflectance = a*R(1,:) + (1-a)*R(5,:);
    else%1,4,11
        transmittance = T(5,:);
        reflectance = R(5,:);       
    end
elseif SHGC>.6 
    b = (SHGC-.6)/.05;
    if U>4.54
        transmittance = b*T(1,:) + (1-b)*T_bdcd;  
        reflectance = b*R(1,:) + (1-b)*R_bdcd; 
    elseif U>3.41
        a = (U-3.41)/1.13;
        transmittance = a*(b*T(1,:) + (1-b)*T_bdcd) + (1-a)*T(5,:);    
        reflectance = a*(R(1,:) + (1-b)*R_bdcd) + (1-a)*R(5,:);
    else
        transmittance = T(5,:);
        reflectance = R(5,:);
    end
elseif SHGC>.55 
    if U>4.54
        transmittance = T_bdcd;  
        reflectance = R_bdcd; 
    elseif U>3.41
        a = (U-3.41)/1.13;
        transmittance = a*T_bdcd + (1-a)*T(5,:);    
        reflectance = a*R_bdcd + (1-a)*R(5,:);
    else
        transmittance = T(5,:);
        reflectance = R(5,:);
    end
elseif SHGC>.5 
    b = (SHGC-.5)/.05;
    if U>4.54
        transmittance = T_bdcd;  
        reflectance = R_bdcd; 
    elseif U>3.41
        a = (U-3.41)/1.13;
        transmittance = a*T_bdcd + (1-a)*(b*T(5,:) + (1-b)*T_fghi);    
        reflectance = a*R_bdcd + (1-a)*(b*R(5,:) + (1-b)*R_fghi);    
    elseif U>1.7
        transmittance = b*T(5,:) + (1-b)*T_fghi;
        reflectance = b*R(5,:) + (1-b)*R_fghi;
    elseif U>1.42
        a = (U-1.42)/0.28;
        transmittance = a*(b*T(5,:) + (1-b)*T_fghi) + (1-a)*T(5,:);    
        reflectance = a*(b*R(5,:) + (1-b)*R_fghi) + (1-a)*R(5,:);   
    else
        transmittance = T(5,:);
        reflectance = R(5,:);
    end
elseif SHGC>.45 
    if U>4.54
        transmittance = T_bdcd;  
        reflectance = R_bdcd; 
    elseif U>3.41
        a = (U-3.41)/1.13;
        transmittance = a*T_bdcd + (1-a)*T_fghi;    
        reflectance = a*R_bdcd + (1-a)*R_fghi;    
    elseif U>1.7
        transmittance = T_fghi;
        reflectance = R_fghi;
    elseif U>1.42
        a = (U-1.42)/0.28;
        transmittance = a*T_fghi + (1-a)*T(5,:);    
        reflectance = a*R_fghi + (1-a)*R(5,:);   
    else
        transmittance = T(5,:);
        reflectance = R(5,:);
    end
elseif SHGC>.35 
    if U>4.54
        b = (SHGC-.3)/.15;
        transmittance = b*T_bdcd + (1-b)*T(4,:);  
        reflectance = b*R_bdcd + (1-b)*R(4,:);  
    elseif U>3.41
        a = (U-3.41)/1.13;
        b = (SHGC-.3)/.15;
        transmittance = a*(b*T_bdcd + (1-b)*T(4,:)) + (1-a)*T_fghi;    
        reflectance = a*(b*R_bdcd + (1-b)*R(4,:)) + (1-a)*R_fghi;    
    elseif U>1.7
        transmittance = T_fghi;
        reflectance = R_fghi;
    elseif U>1.42
        b = (SHGC-.35)/.1;
        a = (U-1.42)/0.28;
        transmittance = a*T_fghi + (1-a)*(b*T(5,:) + (1-b)*T(10,:));    
        reflectance = a*R_fghi + (1-a)*(b*R(5,:) + (1-b)*R(10,:));   
    else
        b = (SHGC-.35)/.1;
        transmittance = b*T(5,:) + (1-b)*T(10,:);
        reflectance = b*R(5,:) + (1-b)*R(10,:);
    end
elseif SHGC>.3 
    if U>4.54
        b = (SHGC-.3)/.15;
        transmittance = b*T_bdcd + (1-b)*T(4,:);  
        reflectance = b*R_bdcd + (1-b)*R(4,:);  
    elseif U>3.41
        a = (U-3.41)/1.13;
        b = (SHGC-.3)/.15;
        transmittance = a*(b*T_bdcd + (1-b)*T(4,:)) + (1-a)*T_fghi;    
        reflectance = a*(b*R_bdcd + (1-b)*R(4,:)) + (1-a)*R_fghi;    
    elseif U>1.7
        transmittance = T_fghi;
        reflectance = R_fghi;
    elseif U>1.42
        a = (U-1.42)/0.28;
        transmittance = a*T_fghi + (1-a)*T(10,:);    
        reflectance = a*R_fghi + (1-a)*R(10,:);   
    else
        transmittance = T(10,:);
        reflectance = R(10,:);
    end    
elseif SHGC>.25 
    b = (SHGC-.25)/.05;
    if U>4.54
        transmittance = T(4,:);  
        reflectance = R(4,:);  
    elseif U>3.41
        a = (U-3.41)/1.13;
        transmittance = a*T(4,:) + (1-a)*(b*T_fghi + (1-b)*T_fh);    
        reflectance = a*R(4,:) + (1-a)*(b*R_fghi + (1-b)*R_fh);    
    elseif U>1.7
        transmittance = b*T_fghi + (1-b)*T_fh;
        reflectance = b*R_fghi + (1-b)*R_fh;
    elseif U>1.42
        a = (U-1.42)/0.28;
        transmittance = a*(T_fghi + (1-b)*T_fh) + (1-a)*T(10,:);    
        reflectance = a*(b*R_fghi + (1-b)*R_fh) + (1-a)*R(10,:);   
    else
        transmittance = T(10,:);
        reflectance = R(10,:);
    end  
else
    if U>4.54
        transmittance = T(4,:);  
        reflectance = R(4,:);  
    elseif U>3.41
        a = (U-3.41)/1.13;
        transmittance = a*T(4,:) + (1-a)*T_fh;    
        reflectance = a*R(4,:) + (1-a)*R_fh;    
    elseif U>1.7
        transmittance = T_fh;
        reflectance = R_fh;
    elseif U>1.42
        a = (U-1.42)/0.28;
        transmittance = a*(T_fh) + (1-a)*T(10,:);    
        reflectance = a*(R_fh) + (1-a)*R(10,:);   
    else
        transmittance = T(10,:);
        reflectance = R(10,:);
    end  
end
end%Ends funtion transmit_reflect_coef