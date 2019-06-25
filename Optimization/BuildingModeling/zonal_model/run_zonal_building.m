function [T_zone,T_surf,m_v_zone,m_sys,Q_transfer,facility,Q_test_plots] = run_zonal_building(building,temperatures,humidity,weather,date,des_day)
%Run full simulation of multi-zone building (with warm up if neccessary)
air_density = 1.204; %kg/m^3 %standard air density (1.204 kg/m3) adjusted for the local barometric pressure (standard barometric pressure corrected for altitude, ASHRAE 1997 HOF pg. 6.1).
if length(weather.DrybulbC)~=(length(date)-1)
    weather = interpolate_weather(weather,date(2:end));
end
n_s = length(date)-1;
n = length(building.zones.name);

%pre-size some logging variables
T_zone = zeros(n_s+1,n);
T_surf = zeros(n_s+1,length(building.ctf.subsurf_names));
m_v_zone = zeros(n_s,n);
Q_transfer = zeros(n_s,n);%energy transfer to zone
m_sys = zeros(n_s,n);%air mass flow to zone
frost = zeros(n_s+1,length(building.cases.name)); %need to make as initial condition that is passed from previous time step
cat = {'heat_elec';'heat_gas';'cool_elec';'fan_elec';'water_gas';'water_elec';'tower_elec';'pumps'};
for i = 1:1:length(cat)
    e_use.(cat{i}) = 0;
    facility.(cat{i}) = zeros(n_s,1);
end
%load given conditions or compute warm-up period
if ~isempty(temperatures)
    T_zone(1,:) = temperatures(1:n)';
    T_surf(1,:) = temperatures(n+1:end)';
    m_v_zone(1,:) = humidity';
    w.zone = humidity;
else
    wu_days = building.site.min_warm_up +.5*(building.site.max_warm_up - building.site.min_warm_up);%warm up days
    dp = nnz(date<(date(1)+1));
    wu_date = date(1:dp+1);%1st 24 hours
    for i = 2:wu_days
        wu_date = [wu_date;linspace(date(2),date(dp+1),dp)';];
    end
    wu_weather = interpolate_weather(weather,wu_date(2:end));
    temperatures = 23*ones(n+length(building.ctf.subsurf_names),1);
    humidity = 0.00825*ones(n,1);%initial humidity of each zone
    [T_z,T_s,m_v] = run_zonal_building(building,temperatures,humidity,wu_weather,wu_date,des_day);
    T_zone(1,:) = T_z(end,:);
    T_surf(1,:) = T_s(end,:);
    m_v_zone(1,:) = m_v(end,:);
    temperatures = [T_z(end,:)';T_s(end,:)';];
    w.zone = m_v_zone(1,:)';
end
[T,air_nodes,plant_nodes] = intial_param(building,temperatures,date(1),des_day);
weather_f = fieldnames(weather);
weather_f = weather_f(~strcmp('filename',weather_f));
t = 1;
t_last = date(t);
while length(date)>t  %using while loop so it can be set as continuous
    dt = (date(t+1) - t_last)*24*3600;
    if dt <0
        dt = dt + 24*3600;%necessary for warm-up period when date jumps back to start of day
    end
    for i = 1:1:length(weather_f)
        weather_now.(weather_f{i}) = weather.(weather_f{i})(t);
    end
    [T,w,air_nodes,plant_nodes,e_use] = reset_loop_param(building,T,w,air_nodes,plant_nodes,e_use,cat,weather_now,date(t+1));
    schedules = load_sched(building.schedule,building.holidays,date(t+1),[],des_day);
    [T,w,central_air,direct_air,mix,plenum,air_nodes,plant_nodes,frost(t+1,:),e_use,gains,window_gain,infiltration] = hvac_loads(building,schedules,air_nodes,T,w,plant_nodes,e_use,frost(t,:),weather_now,date(t+1),dt,des_day);
    %% Alter schedule to increase/decrease temperature & compute alternate loads
    
    %% simulate a time step
    central_flow = central_air.flow;
    direct_flow = direct_air.flow;
    T.direct = direct_air.T;
    T.central = central_air.T;
    w.direct = direct_air.w;
    w.central = central_air.w;
    [T.zone,T.surf,w.zone] = update_model_states_implicit(building,T,gains,window_gain,dt,w,mix,plenum,infiltration',central_flow,direct_flow,weather_now);

    %% update logging variables and move to next time step
    T_zone(t+1,:) = T.zone';
    T_surf(t+1,:) = T.surf';
    m_sys(t,:) = air_density*(central_flow + direct_flow)';
    m_v_zone(t+1,:) = w.zone';
    Q_transfer(t,:) = air_density*(central_flow'.*(1006 + 1860*w.central').*(T.central' - T.zone') + direct_flow'.*(1006 + 1860*w.zone').*(T.direct' - T.zone') + (1006 + 1860*w.zone').*(plenum*T.zone - sum(plenum,2).*T.zone)').*building.zones.multiplier';%Sensible energy transfer to the zone.
    for i = 1:1:length(cat)
        facility.(cat{i})(t) = e_use.(cat{i});
    end
    
    %% calculations just for plotting/validation
    air_density = 1.204; %kg/m^3 %standard air density (1.204 kg/m3) adjusted for the local barometric pressure (standard barometric pressure corrected for altitude, ASHRAE 1997 HOF pg. 6.1).
    cp_zone = (1006 + 1860*w.zone)*air_density;%Specific heat in J/m^3
    cp_air = (1006 + 1860*w.air)*air_density;%Specific heat in J/m^3
    T_surf_interior = T.surf(building.ctf.sur_state1);
    dT = T_surf_interior - building.ctf.map_sur2zone'*T.zone;%interior surface convection = h*A*(Tsur - Tzone) 
    h.interior = natural_convection(building.surfaces.normal,dT,[],building.convection.interior);
    surface2zone = h.interior.*building.surfaces.area.*dT;
    h_windows = window_convection(building.windows,(T.windows_int+273),(T.zone(building.windows.zone_index)+273),cp_zone(building.windows.zone_index));
    window2zone = h_windows.*building.windows.area.*(T.windows_int - T.zone(building.windows.zone_index));

    Q_test_plots.infiltration(t,:) = infiltration.*cp_air.*(T.air_zone - T.zone)';
    Q_test_plots.surface_conv(t,:) = (building.ctf.map_sur2zone*surface2zone + building.ctf.map_win2zone*window2zone)';
    Q_test_plots.net_internal_zone_gain(t,:) =  gains.zone_sensible + gains.zone_latent + Q_test_plots.surface_conv(t,:) + Q_test_plots.infiltration(t,:);
    % sig = 5.67e-8;% W/m^2*K^4 Stephan-Boltzman const
    % n_sur = length(building.surfaces.name);
    % area_i = building.ctf.interior_absorb.*building.surfaces.area_rad*ones(1,n_sur);
    % area_w = building.ctf.interior_absorb.*building.surfaces.area_rad*ones(1,length(building.windows.name));
    % Ti = (T.surf(building.ctf.sur_state1)+273)*ones(1,n_sur);
    % surface2surface = sig*building.zones.view_factors(1:n_sur,1:n_sur).*area_i.*(Ti.^4 - Ti'.^4);%heat transfer from i to j: h_ij*(Ti - Tj) = Q_i -->j
    % surface2window = sig*building.zones.view_factors(1:n_sur,n_sur+1:end).*area_w.*((T_surf_interior+273).^4 - (T.windows_int+273)'.^4);%heat transfer from i to j: h_ij*(Ti - Tj) = Q_i -->j
    % T_surf_exterior = T.surf(building.ctf.sur_state2);
    % [h.exterior,h.sky] = exterior_convection(building,T.surf(building.ctf.sur_state2),T.air_surface,weather,T.sky,building.convection.exterior);
    % q_exterior = (gains.exterior_surface' + h.exterior.*(T.exterior - T_surf_exterior) + h.sky.*(T.sky - T_surf_exterior) + building.ctf.ground_cond.*(T.ground - T_surf_exterior)).*building.ctf.exterior_area;
    % q_interior = gains.interior_surface' + window_gain - surface2zone - sum(surface2surface,2) - sum(surface2window,2);
%     Q_test_plots.interior_surface_heat_transfer(t,:) = surface2zone';
%     Q_test_plots.interior_surface_temperature(t,:) = T_surf_interior';
%     Q_test_plots.exterior_surface_heat_transfer(t,:) = q_exterior';
%     Q_test_plots.exterior_surface_temperature(t,:) = T_surf_exterior';
%     Q_test_plots.window_heat_transfer(t,:) = window2zone';
%     Q_test_plots.window_temperature(t,:) = T.windows_int';
    
    t= t+1;
    t_last = date(t);
end
end%Ends function run_zonal_building

function [T,air_nodes,plant_nodes] = intial_param(building,temperatures,date,des_day)
%initial conditions for some more state variables
n = length(building.zones.name);
T.zone = temperatures(1:n);
T.surf = temperatures(n+1:end);
T.windows = 300*ones(length(building.windows.name),4);
SG = f_index('SimpleGlazingSystem',building.windows.type);
T.windows(SG,3:4) = 0;
T.tank = [];
T.tank_on = false(length(building.water_heater.name),1);
for i = 1:1:length(building.water_heater.name)
    T.tank(i,1) = load_sched(building.schedule,building.holidays,date,building.water_heater.temperature_schedule(i),des_day);
end
if isempty(building.plant_demand_nodes)
    building.plant_demand_nodes.name = [];
end
plant_nodes.demand_temperature = zeros(length(building.plant_demand_nodes.name),1);
plant_nodes.demand_flow = zeros(length(building.plant_demand_nodes.name),1);
plant_nodes.supply_temperature = zeros(length(building.plant_supply_nodes.name),1);
for i = 1:1:length(building.plant_loop.name)
    plant_nodes.demand_temperature(building.plant_demand_nodes.loop == i) = building.plant_loop.exit_temperature(i);
    if strcmpi(building.plant_loop.type{i},'Cooling')
        plant_nodes.supply_temperature(building.plant_supply_nodes.loop == i) = building.plant_loop.exit_temperature(i)+building.plant_loop.temperature_difference(i);
    else%heating loop, water comes back cooler
        plant_nodes.supply_temperature(building.plant_supply_nodes.loop == i) = building.plant_loop.exit_temperature(i)-building.plant_loop.temperature_difference(i);
    end
end
if isempty(building.air_supply_nodes)
    building.air_supply_nodes.name = [];
end
sn = length(building.air_supply_nodes.name);
air_nodes.supply_T_set = zeros(sn,1);
air_nodes.supply_w_set = zeros(sn,1);
air_nodes.supply_flow = zeros(sn,1);
dn = length(building.air_demand_nodes.name);
air_nodes.demand_T_set = zeros(dn,1);
air_nodes.demand_w_set = zeros(dn,1);
air_nodes.demand_flow = zeros(dn,1);
end%Ends function initial_param

function [T,w,air_nodes,plant_nodes,e_use] = reset_loop_param(building,T,w,air_nodes,plant_nodes,e_use,cat,weather,date)
if ~isempty(building.plant_demand_nodes)
    plant_nodes.load =  0*plant_nodes.demand_temperature; %reset to zero, and put in loads from demand side as you get through equipment
    plant_nodes.demand_flow = 0*plant_nodes.demand_flow; %reset to zero, and put in loads from demand side as you get through equipment
    plant_nodes.supply_flow = 0*plant_nodes.supply_temperature;
end
if ~isempty(building.air_supply_nodes)
    air_nodes.supply_T_set = 0*air_nodes.supply_T_set;
    air_nodes.supply_flow = 0*air_nodes.supply_flow;
    air_nodes.demand_T_set = 0*air_nodes.demand_T_set;
    air_nodes.demand_flow = 0*air_nodes.demand_flow;
end
for i = 1:1:length(cat)
    e_use.(cat{i}) = 0;
end
%%pre-compute some ambient conditions
sig = 5.67e-8;% W/m^2*K^4 Stephan-Boltzman const
sky_emissivity = (0.787 + 0.764*log((weather.DewpointC+273)/273)).*(1 + 0.0224*weather.OpqCldtenths - 0.0035*weather.OpqCldtenths.^2 + 0.00028*weather.OpqCldtenths.^3);
ir_intensity = sky_emissivity*sig.*(weather.DrybulbC+273).^4;%W/m^2
T.sky = (ir_intensity/sig).^.25 - 273;
T.ground = building.ground_temperatures(month(date));
w.air = psychometric('P',weather.PressurePa,'Tdp',weather.DewpointC,'w');
T.air_zone = weather.DrybulbC - 0.0065*building.ctf.z_height;%outdoor air temp next to each zone
T.air_surface = weather.DrybulbC - 0.0065*building.ctf.s_height;%outdoor air temp next to each exterior surface
T.exterior = T.air_surface;
T.exterior(f_index('Ground',building.surfaces.exterior.boundary)) = T.ground;
T.surf_old = T.surf;
end%Ends function reset_loop_param