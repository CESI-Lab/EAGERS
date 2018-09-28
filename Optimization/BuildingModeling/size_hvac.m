function building = size_hvac(building,weather,size_info)
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
air_density = 1.225; %kg/m^3

%% demands based on winter design day
dd = 1;
d1 = datenum([2017,building.HVAC.design.month(dd),building.HVAC.design.day(dd)]);
w_date = linspace(d1,d1+1,24/ts+1)';%hourly for a day
w_weather = interpolate_weather(weather,w_date(2:end));
a = building.HVAC.design.drybulb_range(dd)/(max(w_weather.DrybulbC) - min(w_weather.DrybulbC));
w_weather.DrybulbC = a*(w_weather.DrybulbC - min(w_weather.DrybulbC)) + (building.HVAC.design.drybulb_temp(dd) - building.HVAC.design.drybulb_range(dd));
w_weather.DewpointC = ones(24/ts,1)*building.HVAC.design.wetbulb_temp(dd);
w_weather.PressurePa= ones(24/ts,1)*building.HVAC.design.pressure(dd);
w_weather.Wspdms = ones(24/ts,1)*building.HVAC.design.windspeed(dd);
w_weather.Wdirdegrees = ones(24/ts,1)*building.HVAC.design.winddir(dd);
w_weather.OpqCldtenths = ones(24/ts,1)-building.HVAC.design.clearness(dd);

[T_zone,~,~,m_sys_w,Q_sys_w,T_sup_w,T_mix_w,w_mix_w] = zone_simulate(building,[],[],w_weather,w_date,'winter');
zone_flow_w  = roll_window(m_sys_w/air_density,building.HVAC.design.averaging_window/ts);%flow in m^3/s
Q_sys_w_avg  = roll_window(Q_sys_w,building.HVAC.design.averaging_window/ts);
size_param.zone.calc_load(:,1) =  max(0,max(Q_sys_w_avg))';
size_param.zone.user_load(:,1) =  size_param.zone.calc_load(:,1)*building.HVAC.design.heat_sizing_factor;
size_param.zone.calc_flow(:,1) =  max(zone_flow_w)';
size_param.zone.user_flow(:,1) =  size_param.zone.calc_flow(:,1)*building.HVAC.design.heat_sizing_factor;

%% demands based on summer design day
dd = 2;
d2 = datenum([2017,building.HVAC.design.month(dd),building.HVAC.design.day(dd)]);
s_date = linspace(d2,d2+1,24/ts+1)';%hourly for a day
s_weather = interpolate_weather(weather,s_date(2:end));
a = building.HVAC.design.drybulb_range(dd)/(max(s_weather.DrybulbC) - min(s_weather.DrybulbC));
s_weather.DrybulbC = a*(s_weather.DrybulbC - min(s_weather.DrybulbC)) + (building.HVAC.design.drybulb_temp(dd) - building.HVAC.design.drybulb_range(dd));
s_weather.DewpointC = ones(24/ts,1)*building.HVAC.design.wetbulb_temp(dd);
s_weather.PressurePa= ones(24/ts,1)*building.HVAC.design.pressure(dd);
s_weather.Wspdms = ones(24/ts,1)*building.HVAC.design.windspeed(dd);
s_weather.Wdirdegrees = ones(24/ts,1)*building.HVAC.design.winddir(dd);
s_weather.OpqCldtenths = ones(24/ts,1)-building.HVAC.design.clearness(dd);

[T_zone,~,~,m_sys_s,Q_sys_s,T_sup_s,T_mix_s,w_mix_s] = zone_simulate(building,[],[],s_weather,s_date,'summer');
zone_flow_s  = roll_window(m_sys_s/air_density,building.HVAC.design.averaging_window/ts);%flow in m^3/s
Q_sys_s_avg  = roll_window(Q_sys_s,building.HVAC.design.averaging_window/ts);
size_param.zone.calc_load(:,2) =  max(0,max(-Q_sys_s_avg))';
size_param.zone.user_load(:,2) =  size_param.zone.calc_load(:,2)*building.HVAC.design.cooling_sizing_factor;
size_param.zone.calc_flow(:,2) =  max(zone_flow_s)';
size_param.zone.user_flow(:,3) =  size_param.zone.calc_flow(:,2)*building.HVAC.design.cooling_sizing_factor;

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

zone_heat_flow_w = zone_flow_w.*(Q_sys_w_avg>0)*building.HVAC.design.heat_sizing_factor;%flow in m^3/s
zone_cool_flow_w = zone_flow_w.*(Q_sys_w_avg<0)*building.HVAC.design.cooling_sizing_factor;%flow in m^3/s
zone_main_flow_w = zone_heat_flow_w + zone_cool_flow_w;%flow in m^3/s
zone_cool_w = -Q_sys_w_avg.*((-Q_sys_w_avg)>0)*building.HVAC.design.cooling_sizing_factor;
zone_heat_w = Q_sys_w_avg.*(Q_sys_w_avg>0)*building.HVAC.design.heat_sizing_factor;

%%need to add heating to get fresh air up to T_supply?
%%need to remove fan heating load?


zone_heat_flow_s = zone_flow_s.*(Q_sys_s_avg>0)*building.HVAC.design.heat_sizing_factor;%flow in m^3/s
zone_cool_flow_s = zone_flow_s.*(Q_sys_s_avg<0)*building.HVAC.design.cooling_sizing_factor;%flow in m^3/s
zone_main_flow_s = zone_heat_flow_s + zone_cool_flow_s;%flow in m^3/s
zone_cool_s = -Q_sys_s_avg.*((-Q_sys_s_avg)>0)*building.HVAC.design.cooling_sizing_factor;
zone_heat_s = Q_sys_s_avg.*(Q_sys_s_avg>0)*building.HVAC.design.heat_sizing_factor;

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
    loop_cool_w(t,:) = (building.HVAC.zone_2_loop*(zone_cool_w(t,:)'.*building.zones.multiplier))';
    loop_heat_w(t,:) = (building.HVAC.zone_2_loop*(zone_heat_w(t,:)'.*building.zones.multiplier))';
    loop_main_flow_w(t,:) = (building.HVAC.zone_2_loop*((zone_heat_flow_w(t,:) + zone_cool_flow_w(t,:))'.*building.zones.multiplier))';%flow in m^3/s
end


for t = 1:1:24/ts
    loop_cool_s(t,:) = (building.HVAC.zone_2_loop*(zone_cool_s(t,:)'.*building.zones.multiplier))';
    loop_heat_s(t,:) = (building.HVAC.zone_2_loop*(zone_heat_s(t,:)'.*building.zones.multiplier))';
    loop_main_flow_s(t,:) = (building.HVAC.zone_2_loop*((zone_heat_flow_s(t,:) + zone_cool_flow_s(t,:))'.*building.zones.multiplier))';%flow in m^3/s
end
size_param.loop_flow = building.HVAC.zone_2_loop*(size_param.zone_flow.*building.zones.multiplier);


for i =1:1:length(size_param.loop_flow)
    [size_param.loop.user_load(i,1),I] = max(loop_heat_w(:,i));
    size_param.loop.user_flow(i,1) = loop_main_flow_w(I,i);
    size_param.loop.hour_of_day(i,1) = I*ts;
    [size_param.loop.user_load(i,2),I] = max(loop_cool_s(:,i));
    size_param.loop.user_flow(i,2) = loop_main_flow_s(I,i);
    size_param.loop.hour_of_day(i,2) = I*ts;
end

%% ---%% Reset acording to .eio sizing results file
size_param = size_info;
%%  ---- %%%
for i = 1:1:length(building.HVAC.terminals.name)
    if isnan(building.HVAC.terminals.max_flow(i))
        building.HVAC.terminals.max_flow(i) = size_param.zone_flow(building.HVAC.terminals.zone(i));%size if not set initially
    end
end

building.HVAC.design.max_loop_flow = size_param.loop_flow;
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

if any(~strcmp(building.HVAC.terminals.type,'Uncontrolled'))
    building.HVAC.terminals.min_flow = building.HVAC.terminals.max_flow.*building.HVAC.terminals.min_frac;
    reheat_term = nonzeros((1:length(building.HVAC.terminals.name))'.*strcmp(building.HVAC.terminals.type,'Reheat'));
    for j = 1:1:length(reheat_term)
        i = reheat_term(j);%terminal index
        name = building.HVAC.terminals.reheat_coil_name{i};
        k = nonzeros((1:length(building.HVAC.coils.heating.name))'.*strcmpi(name,building.HVAC.coils.heating.name));%Coil index
        al = building.HVAC.terminals.loop(i);%air loop
        switch building.HVAC.coils.heating.type{k}
            case {'Fuel';'Electric'}
                c = nonzeros((1:length(size_info.component.name))'.*strcmpi(name,size_info.component.name));
                building.HVAC.coils.heating.capacity(k) = size_info.component.value(c);
            case 'Water'    
                building.HVAC.coils.heating = heating_coil(building.HVAC.coils.heating,building.HVAC.coils.heating.name{k},'sizing',peak.heat_mix_temp(al),peak.heat_mix_w(al),peak.heat_flow(al)*air_density);
                building.HVAC.terminals.max_flow_water(i) = building.HVAC.coils.heating.rated_water_flow(k);
            otherwise
                disp('need to add missing equipment type in reheat sizing')
        end
    end
end

%%set outdoor air minimum and maximum
for i = 1:1:length(building.HVAC.outdoor_air.name)
    k = nonzeros((1:length(size_param.component.name))'.*strcmpi(building.HVAC.outdoor_air.control(i),size_param.component.name));
    k2 = 1;
    while k2<length(k) && ~strcmpi(size_param.component.field{k(k2)},'Minimum Outdoor Air Flow Rate [m3/s]')
        k2 = k2+1;
    end
    k = k(k2);
    l = building.HVAC.outdoor_air.loop(i);
    if l>0
        building.HVAC.outdoor_air.min_flow(i) = size_param.component.value(k);
        building.HVAC.outdoor_air.max_flow(i) = size_param.loop_flow(l);
    else
        z = building.HVAC.outdoor_air.zone(i);
        building.HVAC.outdoor_air.min_flow(i) = size_param.component.value(k);
        building.HVAC.outdoor_air.max_flow(i) = size_param.zone_flow(z);
    end
end

%loop through equipment on loop
for al = 1:1:length(building.HVAC.zone_2_loop(:,1))%air loop list
    output.cool = peak.cool_load(al);
    output.cool_T = peak.cool_mix_temp(al);
    output.cool_w = peak.cool_mix_w(al);
    output.cool_m_dot = peak.cool_flow(al)*air_density;
    output.cool_des_T = building.HVAC.design.Tsupply_c(al);
    output.cool_des_w = building.HVAC.design.w_supply_c(al);
    output.all_cool = building.HVAC.design.all_outdoor_flow_cooling(al);
    output.heat = peak.heat_load(al);
    output.heat_T = peak.heat_mix_temp(al);
    output.heat_w = peak.heat_mix_w(al);
    output.heat_m_dot = peak.heat_flow(al)*air_density;
    output.nominal_flow = size_param.loop_flow(al);
    
    for j = 1:1:length(building.HVAC.loop.component_name{al})
        type = building.HVAC.loop.component_type{al}{j};
        name = building.HVAC.loop.component_name{al}{j};
        switch type
            case 'AirLoopHVAC:OutdoorAirSystem'
                %currently do nothing, later some version might have damper or economizer
            case 'AirLoopHVAC:UnitaryHeatCool'
                k = nonzeros((1:length(building.HVAC.unitary_sys.name))'.*strcmp(name,building.HVAC.unitary_sys.name));
                building = unitary_sys_size(building,output,k);
            case {'Coil:Cooling:Water:DetailedGeometry';'Coil:Cooling:Water'}
                building.HVAC.coils.cooling = cooling_coil(building.HVAC.coils.cooling,name,'sizing',output.cool_T,output.cool_w,output.cool_m_dot);
            case {'CoilSystem:Cooling:DX';'DX:TwoSpeed';'DX:SingleSpeed'}
                if ~any(strcmpi(name,building.HVAC.coils.cooling.name))
                    name = building.HVAC.coil_system.coil_name{nonzeros((1:length(building.HVAC.coil_system.name))'.*strcmpi(name,building.HVAC.coil_system.name))};
                end
                building.HVAC.coils.cooling = cooling_coil(building.HVAC.coils.cooling,name,'sizing',building.HVAC.curves,output.nominal_flow,output.cool_des_T,output.cool_des_w,output.all_cool);
            case {'Coil:Heating:Fuel';'Coil:Heating:Electric'}
                k = nonzeros((1:length(building.HVAC.coils.heating.name))'.*strcmpi(name,building.HVAC.coils.heating.name));
                if isnan(building.HVAC.coils.heating.capacity(k))
                    building.HVAC.coils.heating.capacity(k) = output.heat;
                end
            case 'Coil:Heating:Water'
                building.HVAC.coils.heating = heating_coil(building.HVAC.coils.heating,name,'sizing',output.heat_T,output.heat_w,output.heat_m_dot);
            case {'Fan:ConstantVolume';'Fan:VariableVolume';'Fan:OnOff';}
                k = nonzeros((1:length(building.HVAC.fans.name))'.*strcmpi(name,building.HVAC.fans.name));
                if isnan(building.HVAC.fans.flow_rate(k))
                    building.HVAC.fans.flow_rate(k) = output.nominal_flow;
                end
            case 'Humidifier:Steam:Electric'
                %does it need sizing?    
            otherwise
                disp('need to add missing equipment type in sizing')
        end     
    end
end
%% size unitary systems not on an air loop
j = length(size_param.loop_flow);
if ~isempty(building.HVAC.unitary_sys)
    for i = 1:1:length(building.HVAC.unitary_sys.name)
        if building.HVAC.unitary_sys.loop(i) == 0
            z = building.HVAC.unitary_sys.zone(i);
            j = j+1;
            output.cool = size_param.zone.user_load(i,2);    
            output.cool_T = interp1(linspace(0,24,24/ts+1),[T_mix_s(1,j);T_mix_s(:,j)],size_param.zone.hour_of_day(i,2));
            output.cool_w = interp1(linspace(0,24,24/ts+1),[w_mix_s(1,j);w_mix_s(:,j)],size_param.zone.hour_of_day(i,2));
            output.cool_m_dot = size_param.zone_flow(z,2)*air_density;
            output.cool_des_T = building.zones.t_supply_c(z);
            output.cool_des_w = building.zones.humid_supply_c(z);
            output.all_cool = false;
            
            output.heat = size_param.zone.user_load(i,1);    
            output.heat_T = interp1(linspace(0,24,24/ts+1),[T_mix_s(1,j);T_mix_s(:,j)],size_param.zone.hour_of_day(i,1));
            output.heat_w = interp1(linspace(0,24,24/ts+1),[w_mix_s(1,j);w_mix_s(:,j)],size_param.zone.hour_of_day(i,1));
            output.heat_m_dot = size_param.zone_flow(z,1)*air_density;
            output.nominal_flow = size_param.zone_flow(z);
            building = unitary_sys_size(building,output,i);
        end
    end
end

end%Ends function size_HVAC

function building = unitary_sys_size(building,output,i)
unit = building.HVAC.unitary_sys;
f = nonzeros((1:length(building.HVAC.fans.name))'.*strcmp(unit.fan_name{i},building.HVAC.fans.name));
if isnan(building.HVAC.fans.flow_rate(f))
    building.HVAC.fans.flow_rate(f) = output.nominal_flow;
end
if strcmp(unit.type{i},'heat_cool')
    name = unit.cool_coil{i};
    switch unit.cool_coil_type{i}
        case {'CoilSystem:Cooling:DX';'Coil:Cooling:DX:TwoSpeed';'Coil:Cooling:DX:SingleSpeed'}
            if ~any(strcmpi(name,building.HVAC.coils.cooling.name))
                name = building.HVAC.coil_system.coil_name{nonzeros((1:length(building.HVAC.coil_system.name))'.*strcmpi(name,building.HVAC.coil_system.name))};
            end
            building.HVAC.coils.cooling = cooling_coil(building.HVAC.coils.cooling,name,'sizing',building.HVAC.curves,output.nominal_flow,output.cool_des_T,output.cool_des_w,output.all_cool);
        case {'Coil:Cooling:Water:DetailedGeometry';'Coil:Cooling:Water'}
            building.HVAC.coils.cooling = cooling_coil(building.HVAC.coils.cooling,name,'sizing',output.cool_T,output.cool_w,output.cool_m_dot);
    end
end
name = unit.heat_coil{i};
switch unit.heat_coil_type{i}
    case {'Coil:Heating:Fuel';'Coil:Heating:Electric'}
        k = nonzeros((1:length(building.HVAC.coils.heating.name))'.*strcmpi(name,building.HVAC.coils.heating.name));
        if isnan(building.HVAC.coils.heating.capacity(k))
            building.HVAC.coils.heating.capacity(k) = output.heat;
        end
    case 'Coil:Heating:Water'
        building.HVAC.coils.heating = heating_coil(building.HVAC.coils.heating,name,'sizing',output.heat_T,output.heat_w,output.heat_m_dot);
end
if strcmp(unit.type{i},'heat_cool') && ~isempty(unit.reheat_name{i})
    disp('need to add dehumid and re-heat capability to unitary systems')
end
building.HVAC.unitary_sys = unit;
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
