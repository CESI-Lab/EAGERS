function building = size_hvac(building,weather,filename)
%% runs but doesn't give same answers as E+ sizing. Use import_sizes instead

%This function sizes HVAC equipment according to strategy laid out by
%reference manual page 501. Process proceeds as:
% 1) simulate the winter and summer design day
% 2) mass flow is calculated from zone load
% 3) outdoor air requirment is calculated
% 4) flows are temporally averaged over a fixed window
% 5) the peak heating/cooling loads and flows are extracted from the design days
% 6) %The user can modify the calculated zone design results by specifying heating and cooling sizing
% factors at the global or zone level or by specifying and actual design heating or cooling zone
% design volumetric flow rate. All of this input is treated as a sizing factor. If the user inputs
% a cooling design volumetric flow rate for a zone it is divided by the calculated cooling design
% volumetric flow rate for the zone to give a zone cooling sizing factor. Note that the user can
% input a zone sizing factor or a zone design flow rate - not both - so there is never a conflict

%flow agregated into air loops for supply side
ts = 1; %timestep in hours for design days
air_density = 1.204; %kg/m^3 %standard air density (1.204 kg/m3) adjusted for the local barometric pressure (standard barometric pressure corrected for altitude, ASHRAE 1997 HOF pg. 6.1).

%% demands based on winter design day
dd = 1;
d1 = datenum([2017,building.hvac.design.month(dd),building.hvac.design.day(dd)]);
w_date = linspace(d1,d1+1,24/ts+1)';%hourly for a day
w_weather = interpolate_weather(weather,w_date(2:end));
a = building.hvac.design.drybulb_range(dd)/(max(w_weather.DrybulbC) - min(w_weather.DrybulbC));
w_weather.DrybulbC = a*(w_weather.DrybulbC - min(w_weather.DrybulbC)) + (building.hvac.design.drybulb_temp(dd) - building.hvac.design.drybulb_range(dd));
w_weather.DewpointC = ones(24/ts,1)*building.hvac.design.wetbulb_temp(dd);
w_weather.PressurePa= ones(24/ts,1)*building.hvac.design.pressure(dd);
w_weather.Wspdms = ones(24/ts,1)*building.hvac.design.windspeed(dd);
w_weather.Wdirdegrees = ones(24/ts,1)*building.hvac.design.winddir(dd);
w_weather.OpqCldtenths = ones(24/ts,1)-building.hvac.design.clearness(dd);

[T_zone,~,~,m_sys_w,Q_sys_w,T_sup_w,T_mix_w,w_mix_w] = zone_simulate(building,[],[],w_weather,w_date,'winter');
zone_flow_w  = roll_window(m_sys_w/air_density,building.hvac.design.averaging_window/ts);%flow in m^3/s
Q_sys_w_avg  = roll_window(Q_sys_w,building.hvac.design.averaging_window/ts);
size_param.zone.calc_load(:,1) =  max(0,max(Q_sys_w_avg))';
size_param.zone.user_load(:,1) =  size_param.zone.calc_load(:,1)*building.hvac.design.heat_sizing_factor;
size_param.zone.calc_flow(:,1) =  max(zone_flow_w)';
size_param.zone.user_flow(:,1) =  size_param.zone.calc_flow(:,1)*building.hvac.design.heat_sizing_factor;

%% demands based on summer design day
dd = 2;
d2 = datenum([2017,building.hvac.design.month(dd),building.hvac.design.day(dd)]);
s_date = linspace(d2,d2+1,24/ts+1)';%hourly for a day
s_weather = interpolate_weather(weather,s_date(2:end));
a = building.hvac.design.drybulb_range(dd)/(max(s_weather.DrybulbC) - min(s_weather.DrybulbC));
s_weather.DrybulbC = a*(s_weather.DrybulbC - min(s_weather.DrybulbC)) + (building.hvac.design.drybulb_temp(dd) - building.hvac.design.drybulb_range(dd));
s_weather.DewpointC = ones(24/ts,1)*building.hvac.design.wetbulb_temp(dd);
s_weather.PressurePa= ones(24/ts,1)*building.hvac.design.pressure(dd);
s_weather.Wspdms = ones(24/ts,1)*building.hvac.design.windspeed(dd);
s_weather.Wdirdegrees = ones(24/ts,1)*building.hvac.design.winddir(dd);
s_weather.OpqCldtenths = ones(24/ts,1)-building.hvac.design.clearness(dd);

[T_zone,~,~,m_sys_s,Q_sys_s,T_sup_s,T_mix_s,w_mix_s] = zone_simulate(building,[],[],s_weather,s_date,'summer');
zone_flow_s  = roll_window(m_sys_s/air_density,building.hvac.design.averaging_window/ts);%flow in m^3/s
Q_sys_s_avg  = roll_window(Q_sys_s,building.hvac.design.averaging_window/ts);
size_param.zone.calc_load(:,2) =  max(0,max(-Q_sys_s_avg))';
size_param.zone.user_load(:,2) =  size_param.zone.calc_load(:,2)*building.hvac.design.cooling_sizing_factor;
size_param.zone.calc_flow(:,2) =  max(zone_flow_s)';
size_param.zone.user_flow(:,3) =  size_param.zone.calc_flow(:,2)*building.hvac.design.cooling_sizing_factor;

%% During sizing the internal radiant and solar gains are essentially removed
%% See page 1710 for a description & 1716 for subtracting and averaging

% If cooling design air flow method is flow/zone, then cooling design air flow rate will be used for
% the design max cooling air flow rate. If cooling design air flow method is design day, then the design
% day calculation will set the design max cooling air flow rate. If cooling design air flow method is
% design day with limit, then the maximum from cooling min flow per area and cooling min flow will
% set a lower limit on the design max cooling air flow rate. In all cases, the maximum from cooling
% min flow per area, cooling min flow, and cooling min flow fraction will set a minimum zone cooling
% air flow rate. In all cases the maximum design cooling air flow rate must be > = to the ventilation
% requirement.

zone_heat_flow_w = zone_flow_w.*(Q_sys_w_avg>0)*building.hvac.design.heat_sizing_factor;%flow in m^3/s
zone_cool_flow_w = zone_flow_w.*(Q_sys_w_avg<0)*building.hvac.design.cooling_sizing_factor;%flow in m^3/s
zone_main_flow_w = zone_heat_flow_w + zone_cool_flow_w;%flow in m^3/s
zone_cool_w = -Q_sys_w_avg.*((-Q_sys_w_avg)>0)*building.hvac.design.cooling_sizing_factor;
zone_heat_w = Q_sys_w_avg.*(Q_sys_w_avg>0)*building.hvac.design.heat_sizing_factor;

%%need to add heating to get fresh air up to T_supply?
%%need to remove fan heating load?


zone_heat_flow_s = zone_flow_s.*(Q_sys_s_avg>0)*building.hvac.design.heat_sizing_factor;%flow in m^3/s
zone_cool_flow_s = zone_flow_s.*(Q_sys_s_avg<0)*building.hvac.design.cooling_sizing_factor;%flow in m^3/s
zone_main_flow_s = zone_heat_flow_s + zone_cool_flow_s;%flow in m^3/s
zone_cool_s = -Q_sys_s_avg.*((-Q_sys_s_avg)>0)*building.hvac.design.cooling_sizing_factor;
zone_heat_s = Q_sys_s_avg.*(Q_sys_s_avg>0)*building.hvac.design.heat_sizing_factor;

size_param.zone_flow = max(max(zone_main_flow_w)',max(zone_main_flow_s)');
for i =1:1:length(size_param.zone_flow)
    [size_param.zone.user_load(i,1),I] = max(zone_heat_w(:,i));
    size_param.zone.user_flow(i,1) = zone_main_flow_w(I,i);
    size_param.zone.hour_of_day(i,1) = I*ts;
    [size_param.zone.user_load(i,2),I] = max(zone_cool_s(:,i));
    size_param.zone.user_flow(i,2) = zone_main_flow_s(I,i);
    size_param.zone.hour_of_day(i,2) = I*ts;
end

%% agregate to loop & use multiplier

for t = 1:1:24/ts
    loop_cool_w(t,:) = (building.hvac.zone_2_loop*(zone_cool_w(t,:)'.*building.zones.multiplier))';
    loop_heat_w(t,:) = (building.hvac.zone_2_loop*(zone_heat_w(t,:)'.*building.zones.multiplier))';
    loop_main_flow_w(t,:) = (building.hvac.zone_2_loop*((zone_heat_flow_w(t,:) + zone_cool_flow_w(t,:))'.*building.zones.multiplier))';%flow in m^3/s
end


for t = 1:1:24/ts
    loop_cool_s(t,:) = (building.hvac.zone_2_loop*(zone_cool_s(t,:)'.*building.zones.multiplier))';
    loop_heat_s(t,:) = (building.hvac.zone_2_loop*(zone_heat_s(t,:)'.*building.zones.multiplier))';
    loop_main_flow_s(t,:) = (building.hvac.zone_2_loop*((zone_heat_flow_s(t,:) + zone_cool_flow_s(t,:))'.*building.zones.multiplier))';%flow in m^3/s
end
size_param.loop_flow = building.hvac.zone_2_loop*(size_param.zone_flow.*building.zones.multiplier);


for i =1:1:length(size_param.loop_flow)
    [size_param.loop.user_load(i,1),I] = max(loop_heat_w(:,i));
    size_param.loop.user_flow(i,1) = loop_main_flow_w(I,i);
    size_param.loop.hour_of_day(i,1) = I*ts;
    [size_param.loop.user_load(i,2),I] = max(loop_cool_s(:,i));
    size_param.loop.user_flow(i,2) = loop_main_flow_s(I,i);
    size_param.loop.hour_of_day(i,2) = I*ts;
end

%% ---%% Reset acording to .eio sizing results file
size_param = parse_eio2(building,filename);
%%  ---- %%%
for i = 1:1:length(building.hvac.terminals.name)
    %%note: in large hotel some AirTerminal:SingleDuct:Uncontrolled are specified with a maximum air flow rate in the .idf, but the air flow rate is different in the sizing .eio file
%     if isnan(building.hvac.terminals.max_flow(i))
        building.hvac.terminals.max_flow(i) = size_param.zone_flow(building.hvac.terminals.zone(i));%size if not set initially
%     end
end

%%need to set coil water temperatures
%coils.cooling.water_inlet_temperature(i,1) = plant.loop.exit_temperature(plant.components.loop(c_num));
% pl = plant.components.loop(nonzeros((1:length(plant.components.name))'.*strcmpi(coils.cooling.name{i},plant.components.name)));%plant loop
% coils.cooling.water_outlet_temperature(i,1) = coils.cooling.water_inlet_temperature(i,1) + plant.loop.temperature_difference(pl);

building.hvac.design.max_loop_flow = size_param.loop_flow;
%% ---%%%%%
for i =1:1:length(size_param.loop_flow)
    peak.cool_load(i) = size_param.loop.user_load(i,2);    
    peak.cool_mix_temp(i) = interp1(linspace(0,24,24/ts+1),[T_mix_s(1,i);T_mix_s(:,i)],size_param.loop.hour_of_day(i,2));
    peak.cool_mix_w(i) = interp1(linspace(0,24,24/ts+1),[w_mix_s(1,i);w_mix_s(:,i)],size_param.loop.hour_of_day(i,2));
    peak.cool_flow(i) = size_param.loop.user_flow(i,2);
    peak.heat_load(i) = size_param.loop.user_load(i,1);
    peak.heat_mix_temp(i) = interp1(linspace(0,24,24/ts+1),[T_mix_w(1,i);T_mix_w(:,i)],size_param.loop.hour_of_day(i,1));
    peak.heat_mix_w(i) = interp1(linspace(0,24,24/ts+1),[w_mix_w(1,i);w_mix_w(:,i)],size_param.loop.hour_of_day(i,1));
    peak.heat_flow(i) = size_param.loop.user_flow(i,1);
end

if any(~strcmp(building.hvac.terminals.type,'Uncontrolled'))
    building.hvac.terminals.min_flow = building.hvac.terminals.max_flow.*building.hvac.terminals.min_frac;
    reheat_term = nonzeros((1:length(building.hvac.terminals.name))'.*strcmp(building.hvac.terminals.type,'Reheat'));
    for j = 1:1:length(reheat_term)
        i = reheat_term(j);%terminal index
        name = building.hvac.terminals.reheat_coil_name{i};
        k = nonzeros((1:length(building.hvac.coils.heating.name))'.*strcmpi(name,building.hvac.coils.heating.name));%Coil index
        al = building.hvac.terminals.loop(i);%air loop
        switch building.hvac.coils.heating.type{k}
            case {'Fuel';'Electric'}
                c = nonzeros((1:length(size_param.component.name))'.*strcmpi(name,size_param.component.name));
                building.coils.heating.capacity(k) = size_param.component.value(c);
            case 'Water'    
                building.coils.heating = heating_coil_sizing(building.coils.heating,building.coils.heating.name{k},peak.heat_mix_temp(al),peak.heat_mix_w(al),peak.heat_flow(al)*air_density);
                building.hvac.terminals.max_flow_water(i) = building.coils.heating.rated_water_flow(k);
            otherwise
                disp('need to add missing equipment type in reheat sizing')
        end
    end
end

%%set outdoor air minimum and maximum
for i = 1:1:length(building.hvac.outdoor_air.name)
    k = nonzeros((1:length(size_param.component.name))'.*strcmpi(building.hvac.outdoor_air.control(i),size_param.component.name));
    if ~isempty(k)
        k2 = 1;
        while k2<length(k) && ~strcmpi(size_param.component.field{k(k2)},'Minimum Outdoor Air Flow Rate [m3/s]')
            k2 = k2+1;
        end
        k = k(k2);
        l = building.hvac.outdoor_air.loop(i);
        if l>0
            building.hvac.outdoor_air.min_flow(i) = size_param.component.value(k);
            building.hvac.outdoor_air.max_flow(i) = size_param.loop_flow(l);
        else
            z = building.hvac.outdoor_air.zone(i);
            building.hvac.outdoor_air.min_flow(i) = size_param.component.value(k);
            building.hvac.outdoor_air.max_flow(i) = size_param.zone_flow(z);
        end
    end
end

%loop through equipment
for j = 1:1:length(building.hvac.components.name)
    type = building.hvac.components.type{j};
    name = building.hvac.components.name{j};
    loop = building.hvac.components.loop(j);
    
    output.cool = peak.cool_load(loop);
    output.cool_T = peak.cool_mix_temp(loop);
    output.cool_w = peak.cool_mix_w(loop);
    output.cool_m_dot = peak.cool_flow(loop)*air_density;
    output.cool_des_T = building.hvac.design.Tsupply_c(loop);
    output.cool_des_w = building.hvac.design.w_supply_c(loop);
    output.all_cool = building.hvac.design.all_outdoor_flow_cooling(loop);
    output.heat = peak.heat_load(loop);
    output.heat_T = peak.heat_mix_temp(loop);
    output.heat_w = peak.heat_mix_w(loop);
    output.heat_m_dot = peak.heat_flow(loop)*air_density;
    output.nominal_flow = size_param.loop_flow(loop);
    
    switch type
        case 'AirLoopHVAC:OutdoorAirSystem'
            %currently do nothing, later some version might have damper or economizer
        case 'AirLoopHVAC:UnitaryHeatCool'
            k = nonzeros((1:length(building.hvac.unitary_sys.name))'.*strcmp(name,building.hvac.unitary_sys.name));
            building = unitary_sys_size(building,output,k);
        case {'Coil:Cooling:Water:DetailedGeometry';'Coil:Cooling:Water'}
            building.coils.cooling = cooling_coil_sizing(building.coils.cooling,name,output.cool_T,output.cool_w,output.cool_m_dot);
        case {'CoilSystem:Cooling:DX';'DX:TwoSpeed';'DX:SingleSpeed'}
            building.coils.cooling = cooling_coil_sizing(building.coils.cooling,name,building.curves,output.nominal_flow,output.cool_des_T,output.cool_des_w,output.all_cool);
        case {'Coil:Heating:Fuel';'Coil:Heating:Electric'}
            k = nonzeros((1:length(building.coils.heating.name))'.*strcmpi(name,building.coils.heating.name));
            if isnan(building.coils.heating.capacity(k))
                building.coils.heating.capacity(k) = output.heat;
            end
        case 'Coil:Heating:Water'
            building.coils.heating = heating_coil_sizing(building.coils.heating,name,output.heat_T,output.heat_w,output.heat_m_dot);
        case {'Fan:ConstantVolume';'Fan:VariableVolume';'Fan:OnOff';}
            k = nonzeros((1:length(building.fans.name))'.*strcmpi(name,building.fans.name));
            if isnan(building.fans.flow_rate(k))
                building.fans.flow_rate(k) = output.nominal_flow;
            end
        case 'Humidifier:Steam:Electric'
            %does it need sizing?    
        otherwise
            disp('need to add missing equipment type in sizing')
    end     
end
%% size unitary systems not on an air loop
j = length(size_param.loop_flow);
if ~isempty(building.hvac.unitary_sys)
    for i = 1:1:length(building.hvac.unitary_sys.name)
        if building.hvac.unitary_sys.loop(i) == 0
            z = building.hvac.unitary_sys.zone(i);
            j = j+1;
            output.cool = size_param.zone.user_load(i,2);    
            output.cool_T = interp1(linspace(0,24,24/ts+1),[T_mix_s(1,j);T_mix_s(:,j)],size_param.zone.hour_of_day(i,2));
            output.cool_w = interp1(linspace(0,24,24/ts+1),[w_mix_s(1,j);w_mix_s(:,j)],size_param.zone.hour_of_day(i,2));
            output.cool_m_dot = size_param.zone.user_flow(z,2)*air_density;
            output.cool_des_T = building.zones.t_supply_c(z);
            output.cool_des_w = building.zones.humid_supply_c(z);
            output.all_cool = false;
            
            output.heat = size_param.zone.user_load(i,1);    
            output.heat_T = interp1(linspace(0,24,24/ts+1),[T_mix_s(1,j);T_mix_s(:,j)],size_param.zone.hour_of_day(i,1));
            output.heat_w = interp1(linspace(0,24,24/ts+1),[w_mix_s(1,j);w_mix_s(:,j)],size_param.zone.hour_of_day(i,1));
            output.heat_m_dot = size_param.zone.user_flow(z,1)*air_density;
            output.nominal_flow = size_param.zone_flow(z);
            building = unitary_sys_size(building,output,i);
        end
    end
end

end%Ends function size_HVAC

function building = unitary_sys_size(building,output,i)
unit = building.hvac.unitary_sys;
f = nonzeros((1:length(building.fans.name))'.*strcmp(unit.fan_name{i},building.fans.name));
if isnan(building.fans.flow_rate(f))
    building.fans.flow_rate(f) = output.nominal_flow;
end
if strcmp(unit.type{i},'heat_cool')
    name = unit.cool_coil{i};
    switch unit.cool_coil_type{i}
        case {'CoilSystem:Cooling:DX';'Coil:Cooling:DX:TwoSpeed';'Coil:Cooling:DX:SingleSpeed'}
            building.coils.cooling = cooling_coil_sizing(building.coils.cooling,name,building.curves,output.nominal_flow,output.cool_des_T,output.cool_des_w,output.all_cool);
        case {'Coil:Cooling:Water:DetailedGeometry';'Coil:Cooling:Water'}
            building.coils.cooling = cooling_coil_sizing(building.coils.cooling,name,output.cool_T,output.cool_w,output.cool_m_dot);
    end
end
name = unit.heat_coil{i};
switch unit.heat_coil_type{i}
    case {'Coil:Heating:Fuel';'Coil:Heating:Electric'}
        k = nonzeros((1:length(building.coils.heating.name))'.*strcmpi(name,building.coils.heating.name));
        if isnan(building.coils.heating.capacity(k))
            building.coils.heating.capacity(k) = output.heat;
        end
    case 'Coil:Heating:Water'
        building.coils.heating = heating_coil_sizing(building.coils.heating,name,output.heat_T,output.heat_w,output.heat_m_dot);
end
if strcmp(unit.type{i},'heat_cool') && ~isempty(unit.reheat_name{i})
    disp('need to add dehumid and re-heat capability to unitary systems')
end
end%Ends function unitary_sys_size

function output = roll_window(input,n)
[n_s,z] = size(input);
a = floor(n/2);
output = zeros(n_s,z);
for t = 1:n_s
    index = (max(1,t-a):min(n_s,t+a-1));
    output(t,:) = mean(input(index,:));
end
end%Ends function roll_window

function size_info = parse_eio2(building,filename)
z_names = building.zones.name;
l_names = building.hvac.loop.name;
multiplier = building.zones.multiplier;
[cat,heading,z_info] = read_eio2(filename);

zone_info = nonzeros((1:length(cat))'.*strcmp('Zone Sizing Information',cat));
z_objects = z_info{zone_info};
sys = nonzeros((1:length(cat))'.*strcmp('System Sizing Information',cat));
s_objects = z_info{sys};
comp = nonzeros((1:length(cat))'.*strcmp('Component Sizing Information',cat));
c_objects = z_info{comp};

zone.flow = zeros(length(z_names),1);
zone.calc_load = zeros(length(z_names),2);
zone.user_load = zeros(length(z_names),2);
zone.calc_flow = zeros(length(z_names),2);
zone.user_flow = zeros(length(z_names),2);
k = nonzeros((1:length(heading{zone_info}))'.*strcmp('User Des Air Flow Rate {m3/s}',heading{zone_info}));
for i = 1:1:length(z_objects(:,1))
    z = nonzeros((1:length(z_names))'.*strcmpi(fix_name(z_objects{i,1}),z_names));
    new_flow = str2double(z_objects{i,k});
    zone.flow(z) = max(zone.flow(z),new_flow/multiplier(z));
    if strcmpi(z_objects{i,2},'heating')
        a = 1;
    else
        a = 2;
    end
    zone.calc_load(z,a) = str2double(z_objects{i,3});
    zone.user_load(z,a) = str2double(z_objects{i,4});
    zone.calc_flow(z,a) = str2double(z_objects{i,5});
    zone.user_flow(z,a) = str2double(z_objects{i,6});
    if isempty(z_objects{i,8})
        zone.hour_of_day(z,a) = 12;
    else
        zone.hour_of_day(z,a) = time_convert(z_objects{i,8});
    end
end

loop_flow = zeros(length(l_names),1);
loop.user_load = zeros(length(l_names),2);
loop.calc_flow = zeros(length(l_names),2);
loop.user_flow = zeros(length(l_names),2);
k = nonzeros((1:length(heading{sys}))'.*strcmp('User Des Air Flow Rate [m3/s]',heading{sys}));
for i = 1:1:length(s_objects(:,1))
    l = nonzeros((1:length(l_names))'.*strcmpi(s_objects{i,1},l_names));
    new_flow = str2double(s_objects{i,k});
    loop_flow(l) = max(loop_flow(l),new_flow);
    if strcmpi(s_objects{i,2},'heating')
        a = 1;
    elseif strcmpi(s_objects{i,2},'cooling')
        a = 2;
    end
    loop.user_load(l,a) = str2double(s_objects{i,4});
    loop.calc_flow(l,a) = str2double(s_objects{i,5});
    loop.user_flow(l,a) = str2double(s_objects{i,6});
    loop.hour_of_day(l,a) = time_convert(s_objects{i,8});
end

for i = 1:1:length(c_objects(:,1))
    component.type(i,1) = c_objects(i,1);
    component.name(i,1) = c_objects(i,2);
    component.field(i,1) = c_objects(i,3);
    component.value(i,1) = str2double(c_objects{i,4});
end
size_info.zone = zone; size_info.loop = loop; size_info.component = component;
end%ends function parse_eio2

function [cat,heading,z_info] = read_eio2(filename)
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
        cat_i = nonzeros((1:length(cat))'.*strcmp(rmv_spaces(lines{i}(1:com(1)-1)),cat));
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
end%Ends function read_eio2


function hod = time_convert(str)
col = strfind(str,':');
h = str2double(str(col(1)-3:col(1)-1));
m = str2double(str(col(1)+1:col(1)+2));
s = str2double(str(col(2)+1:end));
hod = h+m/60+s/3600;
end%Ends function time_convert

function name = rmv_spaces(name)
while ~isempty(name) && strcmp(name(1),' ')
    name = name(2:end);
end
while ~isempty(name) && strcmp(name(end),' ')
    name = name(1:end-1);
end
end%Ends function rmv_spaces

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