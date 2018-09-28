function varargout = zone_simulate(building,temperatures,humidity,weather,date,des_day)
%Energy balance for multi-zone building
%h refers to convection coefficients.
if length(weather.DrybulbC)~=(length(date)-1)
    weather = interpolate_weather(weather,date(2:end));
end
dt = (date(2:end) - date(1:end-1))*3600*24;
dt(dt<=0) = min(dt(dt>0));%allows date to be the same day multiple times for sizing
n_s = length(dt);

P = weather.PressurePa/1000; % atmospheric pressure (kPa)
T_dp = 273 + weather.DewpointC;
P_H2O_dp = exp((-5.8002206e3)./T_dp + 1.3914993 - 4.8640239e-2*T_dp + 4.1764768e-5*T_dp.^2 - 1.4452093e-8*T_dp.^3 + 6.5459673*log(T_dp))/1000; %saturated water vapor pressure at dewpoint
m_v_air = .621945*(P_H2O_dp./(P-P_H2O_dp));%mass fraction of water in air 
spec_heat_air = 1006 + 1860*m_v_air; %J/kg*K

n = length(building.zones.name);
T_zone = zeros(n_s+1,n);
T_surf = zeros(n_s+1,length(building.ctf.subsurf_names));
m_v_zone = zeros(n_s,n);
T_windows_int = 23*ones(n_s,length(building.windows.name));
T_windows_ext = 23*ones(n_s,length(building.windows.name));
Q_windows = zeros(n_s,n);
Q_sys = zeros(n_s,n);%used for sizing
Q_transfer = zeros(n_s,n);%energy transfer to zone
T_supply = zeros(n_s,n);%Supply air temperature to zone
m_sys = zeros(n_s,n);%air mass flow to zone
T_mixed_air = zeros(n_s,length(building.HVAC.loop.name)+nnz(building.HVAC.unitary_sys.zone));
w_mixed_air = zeros(n_s,length(building.HVAC.loop.name)+nnz(building.HVAC.unitary_sys.zone));
heat_elec = zeros(n_s,1);
heat_gas =  zeros(n_s,1);
cool_elec = zeros(n_s,1);
fans = zeros(n_s,length(building.HVAC.fans.name));
zone_exhaust = zeros(n_s,n);
fresh_air_flow = zeros(n_s,n);
if ~isempty(temperatures)
    T_zone(1,:) = temperatures(1:n);
    T_surf(1,:) = temperatures(n+1:end);
    m_v_zone(1,:) = humidity';
else
    wu_days = 4;%warm up days
    dp = nnz(date<(date(1)+1));
    wu_date = date(1:dp+1);%1st 24 hours
    for i = 2:wu_days
        wu_date = [wu_date;linspace(date(2),date(dp+1),dp)';];
    end
    wu_weather = interpolate_weather(weather,wu_date(2:end));
    temperatures = 23*ones(n+length(building.ctf.subsurf_names),1);
    humidity = 0.00825*ones(n,1);%initial humidity of each zone
    [T_z,T_s,m_v,m_dot] = zone_simulate(building,temperatures,humidity,wu_weather,wu_date,des_day);
    T_zone(1,:) = T_z(end,:);
    T_surf(1,:) = T_s(end,:);
    m_v_zone(1,:) = m_v(end,:);
    m_sys(1,:) = m_dot(end,:);
end
occupancy = zeros(n_s,n);
mixing = zeros(n,n,n_s);
infiltration = zeros(n_s,n);
Q_zone_gain = zeros(n_s,n);
W_zone_gain = zeros(n_s,n);
Q_surf_absorb = zeros(n_s,length(building.ctf.capacitance));%interior, intermediate and exterior surfaces

air_density = 1.225; %kg/m^3
sig = 5.67e-8;% W/m^2*K^4 Stephan-Boltzman const
Sky_emissivity = (0.787 + 0.764*log((weather.DewpointC+273)/273)).*(1 + 0.0224*weather.OpqCldtenths - 0.0035*weather.OpqCldtenths.^2 + 0.00028*weather.OpqCldtenths.^3);
ir_intensity = Sky_emissivity*sig.*(weather.DrybulbC+273).^4;%W/m^2
T_sky = (ir_intensity/sig).^.25 - 273;
[~,mon,~] = datevec(date);
T_ground = building.ground_temperatures(mon);
cos_phi_vert = building.surfaces.normal(:,3)./(building.surfaces.normal(:,1).^2 + building.surfaces.normal(:,2).^2 + building.surfaces.normal(:,3).^2).^.5; %portion of surface pointed upwards

tolerance = 1e-2;
for t = 1:1:n_s
    prof_now = load_sched(building.schedule,building.holidays,date(t+1),[],des_day);
    Tair_z = weather.DrybulbC(t) - 0.0065*building.ctf.z_height;%outdoor air temp next to each zone
    Tair_s = weather.DrybulbC(t) - 0.0065*building.ctf.s_height;%outdoor air temp next to each exterior surface
    [~,~,occupancy(t,:),mixing(:,:,t),infiltration(t,:),Q_zone_gain(t,:),Q_surf_absorb(t,:),...
        W_zone_gain(t,:),internal_irrad_window,T_set_now,cos_phi_windows] = ...
        zone_loads(building,date(t+1),prof_now,weather.DNIWm2(t),weather.DHIWm2(t),T_zone(t,:),m_v_zone(t,:),dt(t),Tair_z');

%%%Part A, iterate to find T_supply and M_supply that achieves temperature constraints
    T_zone_est = T_zone(t,:)';
    T_surf_est = T_surf(t,:)';
    m_v_zone_est = m_v_zone(t,:)';
    zone_est_error = ones(1,length(T_zone_est));
    count = 0;
    while count<3 && any(abs(zone_est_error)>tolerance)%replace with tolerence on T_zone_est        
        net_mixing = plenum_mixing(m_sys(t,:)',mixing(:,:,t),infiltration(t,:)',building.HVAC,building.zones.name);
        [T_windows_int(t,:),T_windows_ext(t,:),Q_windows(t,:),h_windows] = ...
            window_temperature(building,cos_phi_windows',T_zone_est+273,T_surf_est(building.ctf.sur_state1)+273,Tair_z+273,T_sky(t)+273,spec_heat_air(t),weather.DNIWm2(t),weather.DHIWm2(t),weather.Wdirdegrees(t),weather.Wspdms(t),internal_irrad_window,T_windows_int(t,:)'+273,T_windows_ext(t,:)'+273);
        h_rad = internal_radiation(building.zones.view_factors,[building.surfaces.area_rad;building.windows.area],[T_surf_est(building.ctf.sur_state1);T_windows_int(t,:)']+273);
        h = natural_convection(building.surfaces.normal,T_surf_est(building.ctf.sur_state1) - building.ctf.map_sur2zone'*T_zone_est,cos_phi_vert);%interior surface convection = h*A*(Tsur - Tzone) 
        h(strcmp(building.surfaces.type,'InternalMass')) = 0;
        [h_s,h_sky,T_exterior] = exterior_convection(building.surfaces,building.ctf.sur_state2,T_surf_est,Tair_s,weather.Wdirdegrees(t),weather.Wspdms(t),T_sky(t),T_ground(t),cos_phi_vert);%exterior surface convection terms
        
        Q_mixing = air_density*(1006 + 1860*m_v_zone_est).*(net_mixing*T_zone_est - sum(net_mixing,2).*T_zone_est);%heat into zone from another zone
        Q_infiltration = air_density*infiltration(t,:)'.*spec_heat_air(t).*(Tair_z - T_zone_est);%heat into zone from air infiltration 
        Q_surfaces = -building.ctf.map_sur2zone*(h.*(building.ctf.map_sur2zone'*T_zone_est - T_surf_est(building.ctf.sur_state1)).*building.surfaces.area);%heat transfer to zone from walls/floor/ceiling
        [Q_zone_Tcool,Q_zone_Theat] = ...
            zone_balance(Q_zone_gain(t,:)',Q_windows(t,:)',Q_mixing,Q_infiltration,Q_surfaces,building.surfaces.area,building.ctf,T_surf_est,T_set_now,T_zone(t,:)',m_v_zone_est,h,building.zones.volume,dt(t));%%predict zone energy exchange to achieve T_set

        [Q_sys(t,:),m_sys(t,:),T_supply(t,:),T_loop_supply,m_v_supply,fresh_air_zone,zone_exhaust(t,:),mixed_air,non_loop_mixed_air,Q_terminal] = ...
            zone_ideal_HVAC(building,prof_now,net_mixing,infiltration(t,:)',Q_zone_Tcool,Q_zone_Theat,T_set_now,T_zone_est,m_v_zone_est,weather.DrybulbC(t),m_v_air(t),occupancy(t,:)',des_day);
        
        [T_z_new,T_surf_est,m_v_zone_est] = zone_surface_temperature(building,T_zone(t,:)',T_surf(t,:)',...
            Q_surf_absorb(t,:)',h_rad(:,1:length(h)),h,h_s,h_sky,T_exterior,T_sky(t),T_ground(t),...
            dt(t),m_v_zone(t,:)',net_mixing,infiltration(t,:)',spec_heat_air(t),m_sys(t,:)',Q_zone_gain(t,:)',...
            Tair_z,T_windows_int(t,:)',T_supply(t,:)',m_v_supply,h_windows,building.zones.volume,W_zone_gain(t,:)',m_v_air(t));
        zone_est_error = T_z_new - T_zone_est;
        T_zone_est = T_z_new;
        count = count+1;
    end
    T_mixed_air(t,:) = [mixed_air.T;non_loop_mixed_air.T]';
    w_mixed_air(t,:) = [mixed_air.w;non_loop_mixed_air.w]';
    fresh_air_flow(t,:) = fresh_air_zone';
%%%Part B, compute with real HVAC system
    %%%----     %
    if isempty(des_day)
        HVAC_out = manage_air_loops(building,prof_now,m_sys(t,:)',T_supply(t,:)',T_loop_supply,mixed_air,non_loop_mixed_air,Q_terminal,weather.DrybulbC(t));
        T_supply(t,:) = HVAC_out{1}; m_v_supply = HVAC_out{2}; heat_elec(t,:) = HVAC_out{3}; heat_gas(t,:) = HVAC_out{4}; cool_elec(t,:) = HVAC_out{5}; fans(t,:) = HVAC_out{6};
        
        [a,b,c] = zone_surface_temperature(building,T_zone(t,:)',T_surf(t,:)',...
            Q_surf_absorb(t,:)',h_rad(:,1:length(h)),h,h_s,h_sky,T_exterior,T_sky(t),T_ground(t),...
            dt(t),m_v_zone(t,:)',net_mixing,infiltration(t,:)',spec_heat_air(t),m_sys(t,:)',Q_zone_gain(t,:)',...
            Tair_z,T_windows_int(t,:)',T_supply(t,:)',m_v_supply,h_windows,building.zones.volume,W_zone_gain(t,:)',m_v_air(t));
        T_zone(t+1,:) = a';
        T_surf(t+1,:) = b';
        m_v_zone(t+1,:) = c';
        Q_transfer(t,:) = m_sys(t,:).*(1006 + 1860*m_v_supply').*(T_supply(t,:) - T_zone(t+1,:));%Sensible energy transfer to the zone. 
    else
        T_zone(t+1,:) = T_zone_est';
        T_surf(t+1,:) = T_surf_est';
        m_v_zone(t+1,:) = m_v_zone_est';
    end
    %%%%____%%
    
    %%need fix so that when raining outside exterior surface temperature = T wet bulb (page 91 of manual)

    Q_test_plots.Infiltration(t,:) = air_density*infiltration(t,:).*spec_heat_air(t).*(Tair_z' - T_zone(t+1,:));
    Q_test_plots.zone_est_error(t,:) = zone_est_error;
    Q_test_plots.surface_conv(t,:) = building.ctf.map_sur2zone*(h.*building.surfaces.area.*(T_surf(t+1,building.ctf.sur_state1)' - building.ctf.map_sur2zone'*T_zone(t+1,:)'));
    Q_test_plots.air_balance_convective(t,:) = Q_zone_gain(t,:);
    if t<n_s
        m_sys(t+1,:) = m_sys(t,:);%initial guess for next time step
    end
end
if isempty(des_day)
    varargout = {T_zone,T_surf,m_v_zone,m_sys,Q_transfer,heat_elec,heat_gas,cool_elec,Q_test_plots,T_windows_int};
else
    varargout = {T_zone,T_surf,m_v_zone,m_sys,Q_sys,T_supply,T_mixed_air,w_mixed_air};
end
end%Ends function zone_simulate

function [Q_zone_Tset_cool,Q_zone_Tset_heat,Q_surfaces,Q_infiltration,Q_mixing] = ...
    zone_balance(Q_internal,Q_windows,Q_mixing,Q_infiltration,Q_surfaces,area,ctf,T_surf,T_set,T_zone,m_v_zone,h,volume,dt)
%% Predicts the zone energy exchange that achieves the temperature setpoints %% The actual HVAC cooling/heating is different due to recirculation & fresh air
air_density = 1.225; %kg/m^3
spec_heat_zone = 1006 + 1860*m_v_zone; %J/kg*K
%%assume T_zone is at cooling setpoint    
Q_transient = (T_zone - T_set.cool').*air_density.*spec_heat_zone.*volume./dt;%zone heat capacitance*deltaT/dt = power (W)
Q_sur = h.*(ctf.map_sur2zone'*T_set.cool' - T_surf(ctf.sur_state1)).*area;%heat transfer rates at each surface with an interior zone in W (positive is heat into surface)
Q_surfaces_cool = -ctf.map_sur2zone*Q_sur;%heat into the zone from the surface heat transfer
Q_zone_Tset_cool = -Q_internal - Q_windows - Q_surfaces_cool - Q_infiltration - Q_mixing - Q_transient;%Energy provided to the zone in W to achieve T_set.cool, + is heating
Q_zone_Tset_cool(T_set.no_hvac) = 0;

%%assume T_zone is at heating setpoint
Q_transient = (T_zone - T_set.heat').*air_density.*spec_heat_zone.*volume./dt;%zone heat capacitance*deltaT/dt = power (W)
Q_sur = h.*(ctf.map_sur2zone'*T_set.heat' - T_surf(ctf.sur_state1)).*area;%heat transfer rates at each surface with an interior zone in W (positive is heat into surface)
Q_surfaces_heat = -ctf.map_sur2zone*Q_sur;%heat into the zone from the surface heat transfer
Q_zone_Tset_heat = -Q_internal - Q_windows - Q_surfaces_heat - Q_infiltration - Q_mixing - Q_transient;%Energy provided to the zone in W to achieve T_set.heat, + is heating
Q_zone_Tset_heat(T_set.no_hvac) = 0;

%% need to add prediction of humidification/dehumidification load
end%Ends function zone_balance

function [T_int,T_ext,Q_windows,h_i] = window_temperature(building,cos_phi_sun,T_zone,T_surf,T_air_zone,T_sky,cp_air,solar_normal,solar_diffuse,w_dir,w_speed,internal_irrad,T_int,T_ext)
%% all temperatures in Kelvin for this function
windows = building.windows;
%% transmit/reflect for simple glazing
% %option 1
% transmit = windows.solar_transmittance.*(windows.transmittance(:,1).*cos_phi.^4 + windows.transmittance(:,2).*cos_phi.^3 + windows.transmittance(:,3).*cos_phi.^2 + windows.transmittance(:,4).*cos_phi + windows.transmittance(:,5));
% reflect = windows.solar_reflectance.*(windows.reflectance(:,1).*cos_phi.^4 + windows.reflectance(:,2).*cos_phi.^3 + windows.reflectance(:,3).*cos_phi.^2 + windows.reflectance(:,4).*cos_phi + windows.reflectance(:,5));
%option 2 page 295
transmit = windows.solar_transmittance.*cos_phi_sun.*(1+(0.768+0.817*windows.solar_heat_gain.^4).*(1-cos_phi_sun.^2).^(3/2));
f1 = (((2.403*cos_phi_sun - 6.192).*cos_phi_sun + 5.625).*cos_phi_sun - 2.095).*cos_phi_sun + 1;
f2 = (((-1.188*cos_phi_sun + 2.022).*cos_phi_sun + 0.137).*cos_phi_sun - 1.720).*cos_phi_sun;
reflect = windows.solar_reflectance.*(f1 + f2.*windows.solar_heat_gain.^.5)./(0.7413 - (0.7396*windows.solar_heat_gain.^.5));
%% need transmit/reflec as fcn agle for regular glazing


windows.visible_absorptance = max(0,1 - transmit - reflect);%%solve equation on page 290 for absorptance, for single layer it is what is not transmitted or reflected
diffuse_absorptance = windows.visible_absorptance;%cant find any mention describing how to get this paramater for simple window model
T_new = T_int*zeros(1,2);
air_density = 1.225; %kg/m^3
g = 9.81; %m/s^2 gravity
sig = 5.67e-8;% W/m^2*K^4 Stephan-Boltzman const

T_o = T_air_zone(windows.zone_index);
T_i = T_zone(windows.zone_index);
S1 = .5*(solar_normal*cos_phi_sun.*windows.visible_absorptance + solar_diffuse*diffuse_absorptance + internal_irrad.visible'.*windows.visible_absorptance);
S2 = S1 + windows.thermal_absorptance.*internal_irrad.long_wave';

cos_phi_vert = building.windows.normal(:,3)./(building.windows.normal(:,1).^2 + building.windows.normal(:,2).^2 + building.windows.normal(:,3).^2).^.5; %portion of window pointed upwards
F_sky = 0.5*(1+cos_phi_vert);
leward = (cos(w_dir)*windows.normal(:,1) +  sin(w_dir)*windows.normal(:,1))<0;%oreintation relative to wind (surface normal is outward facing)
error = 1;
count = 0;
while error>1e-1 && count<10
    Eo = windows.thermal_absorptance.*windows.emittance_front.*sig.*((0.5*(1-cos_phi_vert)).*T_o.^4 ...
        + F_sky.^1.5.*T_sky.^4 + F_sky.*(1-F_sky.^.5).*T_o.^4);
    Ei = windows.thermal_absorptance.*windows.emittance_back.*sig.*(building.zones.view_factors(length(T_surf)+1:end,:)*[T_surf.^4;T_int.^4]);
    h_n_ext = natural_convection(windows.normal, T_ext - T_o,cos_phi_vert);
    h_o = (h_n_ext.^2 + (3.26*w_speed^.89)^2).^.5;
    h_o(leward) = (h_n_ext(leward).^2 + (3.55*w_speed^.617)^2).^.5;

    deltaT_interior = (T_int - T_i);
    T_mf = T_i + 0.25*(T_int - T_i);
    mu = 3.723e-6 + (4.94e-8)*T_mf;
    lambda = 2.873e-3 + (7.76e-5)*T_mf;
    rayleigh = air_density^2*windows.height.^3*g*cp_air.*abs(deltaT_interior)./(T_mf.*mu.*lambda);
    %need to add correlations for changing tilt angle
    nusselt = 0.13*(rayleigh).^(1/3);%cold window 
    nusselt(deltaT_interior>0) = 0.58*(min(rayleigh(deltaT_interior>0),1e11)).^(1/5);%hot window 
    h_i = nusselt.*(lambda./windows.height);
    
    h_r1 = windows.emittance_front*sig.*T_ext.^3;
    h_r2 = windows.emittance_back*sig.*T_int.^3;
    k1 = windows.thermal_conductivity;
    for i = 1:1:length(windows.name)
        A = [-h_r1(i)-k1(i)-h_o(i), k1(i); k1(i), -h_r2(i)-k1(i)-h_i(i)];
        B = [-Eo(i) - S1(i) - h_o(i)*T_o(i); -Ei(i) - S2(i) - h_i(i)*T_i(i);];
        T_new(i,:) = (A\B)';
    end
    error = max(max(abs(T_new(:,1) - T_ext)),max(abs(T_new(:,2) - T_int)));
    T_ext = (T_new(:,1) + T_ext)/2;
    T_int = (T_new(:,2) + T_int)/2;
    count = count+1;
end
Q_windows = (building.ctf.map_win2zone*(h_i.*(T_int - T_i).*windows.area))';
T_int = T_int - 273;%convert back to Celcius 
T_ext = T_ext - 273;
end%Ends function window_temperature

function mixing = plenum_mixing(m_dot,mixing,infil,hvac,zone_names)
zone_inlet_flow = sum(mixing,2) - sum(mixing,1)' + infil + m_dot;
for i = 1:1:length(hvac.plenum.name)
    z = nonzeros((1:length(zone_names))'.*strcmpi(hvac.plenum.zone{i},zone_names));
    hvac.zone_2_plenum(i,i) = 0; %don't count m_sys or infiltration into a plenum zone as mixing from another zone
    mixing(z,:) = mixing(z,:) + hvac.zone_2_plenum(i,:).*zone_inlet_flow';%flow from column j to row i
end
end%Ends function plenum_mixing

function h = natural_convection(normal,dT,cos_phi)
%%TARP method
vertical = abs(normal(:,3))<1e-3; %vertical surface
upwards = (dT<0 & normal(:,3)>0) | (dT>0 & normal(:,3)<0);
downwards = ~vertical & ~upwards;
h = 1.31*(max(0.01,abs(dT)).^(1/3));%avoid dT = 0
h_up = 9.482/1.31*h./(7.283-abs(cos_phi));
h_down = 1.81/1.31*h./(1.382+abs(cos_phi));
h(upwards) = h_up(upwards);
h(downwards) = h_down(downwards);
end%ends function natural_convection

function [h_s,h_sky,T_exterior] = exterior_convection(surface,sur_state2,Tsurf,Tair,w_dir,w_speed,T_sky,T_ground,cos_phi)
sig = 5.67e-8;% W/m^2*K^4 Stephan-Boltzman const
n_sur = length(surface.name);
ext = nonzeros((1:n_sur)'.*(strcmp(surface.boundary,'Outdoors') | strcmp(surface.boundary,'Ground')));
n_ext = length(ext);
h_s = zeros(n_ext,1);
T_exterior = zeros(n_ext,1);
T_exterior_sur = Tsurf(sur_state2);

F_grnd = 0.5*(1-cos_phi(ext));
F_sky = 0.5*(1+cos_phi(ext));
beta = F_sky.^.5;
h_n = natural_convection(surface.normal((ext),:),T_exterior_sur - Tair,cos_phi(ext));
h_glass = (h_n.^2 + (3.26*w_speed^.89)^2).^.5;
leward = (cos(w_dir)*surface.normal(ext,1) +  sin(w_dir)*surface.normal(ext,1))<0;%oreintation relative to wind (surface normal is outward facing)
h_glass(leward) = (h_n(leward).^2 + (3.55*w_speed^.617)^2).^.5;
h_wind = h_n + surface.roughness_factor(ext).*(h_glass-h_n);

h_ground = surface.absorptance.thermal(ext,2)*sig.*F_grnd.*((T_exterior_sur+273).^4 - (Tair+273).^4)./(T_exterior_sur - Tair);
h_ground(T_exterior_sur == Tair) = 4*Tair(T_exterior_sur == Tair).^3;%avoid nan by using derivative
h_sky = surface.absorptance.thermal(ext,2)*sig.*F_sky.*beta.*((T_exterior_sur+273).^4 - (T_sky+273).^4)./(T_exterior_sur - T_sky);
h_sky(T_exterior_sur == T_sky) = 4*Tair(T_exterior_sur == T_sky).^3;%avoid nan by using derivative
h_air = surface.absorptance.thermal(ext,2)*sig.*F_sky.*(1-beta).*((T_exterior_sur+273).^4 - (Tair+273).^4)./(T_exterior_sur - Tair);
h_air(T_exterior_sur == Tair) = 4*Tair(T_exterior_sur == Tair).^3;%avoid nan by using derivative

a = strcmp(surface.boundary(ext),'Outdoors');
h_s(a) = h_wind(a) + h_air(a) + h_ground(a);
T_exterior(a) = Tair(a);
b = strcmp(surface.boundary(ext),'Ground');
T_exterior(b) = T_ground;
end%Ends function exterior convection

function h_rad = internal_radiation(vf,area,Temp)
%Computes linearized radiation heat transfer coefficients for long-wave radiation between internal surfaces
sig = 5.67e-8;% W/m^2*K^4 Stephan-Boltzman const
n = length(Temp);
area_i = area*ones(1,n);
Ti = Temp*ones(1,n);
h_rad = sig*vf.*area_i.*(Ti.^4 - Ti'.^4)./(Ti-Ti');%heat transfer coefficient from i to j: h_ij*(Ti - Tj) = Q_i -->j
alt = 4*sig*area_i.*vf.*Ti.^3;%avoid nan by using derivative
h_rad(Ti==Ti') = alt(Ti==Ti');
end%Ends function internal_radiation

%% alternative to seperately calculating surface and zone temperatures
function [T_zone,T_surf,m_v_zone] = zone_surface_temperature(building,T_zone,T_surf,Q_surf_absorb,h_rad,h,h_s,h_sky,T_exterior,T_sky,T_ground,dt,m_v_zone,mixing,infiltration,spec_heat_air,m_sys,Q_internal,Tair_z,T_windows,T_supply,m_v_supply,h_windows,volume,W_zone_gain,m_v_air)
%% Analytic equation soln (eqn 2.11 in engineering reference)
% grouping all temperature and humidity states as (T_zone;T_surf;m_v_zone]
z = length(T_zone);
s = length(T_surf);
spec_heat_supply = 1006 + 1860*m_v_supply; %J/kg*K
air_density = 1.225; %kg/m^3
spec_heat_zone = 1006 + 1860*m_v_zone; %J/kg*K
C_T = 1*building.zones.volume;%sensible heat multiplier

%%m*Cp & h*A*T terms for analytic solution  
m_Cp = air_density*(spec_heat_zone.*sum(mixing,2) + infiltration.*spec_heat_air) + m_sys.*spec_heat_supply;
h_A_window = building.ctf.map_win2zone*(h_windows.*building.windows.area);
h_A = building.ctf.map_sur2zone*(h.*building.surfaces.area) + h_A_window;
den_z = h_A + m_Cp;% convection + (inter-zone mixing + infiltration) + HVAC
den_s = sum(-building.ctf.R_neg,2) + building.ctf.interior_surface*(h.*building.surfaces.area + sum(h_rad,1)') + building.ctf.exterior_surface*((h_s + h_sky).*building.ctf.exterior_area + building.ctf.ground_cond);
den_w = air_density*(sum(mixing,2) + infiltration) + m_sys;
den_w(den_w==0) = min(den_w(den_w~=0));%adjustment for zones with zero air flow of any kind into them

%combine
Cap = [air_density*spec_heat_zone.*C_T;building.ctf.capacitance;air_density*volume];%capacitance [zone_thermal;surface;zone_moisture];
den = [den_z;den_s;den_w];
sub_steps = max([10,2*building.surfaces.i_nodes,round(dt/360)]);
for i = 1:1:sub_steps %loop through the wall at a faster time freq, so that internal temperatures distribute
    m_Cp_T = air_density*(spec_heat_zone.*(mixing*T_zone) + infiltration.*spec_heat_air.*Tair_z) +  m_sys.*spec_heat_supply.*T_supply;%flows in m^3/s converted to kg/s*J/(kg*K)*K = W 
    h_A_T = building.ctf.map_sur2zone*(h.*building.surfaces.area.*T_surf(building.ctf.sur_state1)) + building.ctf.map_win2zone*(h_windows.*building.windows.area.*T_windows);
    term1_z = Q_internal + h_A_T + m_Cp_T;% internal gains + convection*T_surface + (inter-zone mixing + infiltration)*T_source + HVAC*T_supply
    term1_s = Q_surf_absorb + building.ctf.interior_surface*(h.*(building.ctf.map_sur2zone'*T_zone).*building.surfaces.area + h_rad'*[T_surf(building.ctf.sur_state1);T_windows]) + building.ctf.R_pos*T_surf + building.ctf.exterior_surface*((h_s.*T_exterior + h_sky.*T_sky).*building.ctf.exterior_area + building.ctf.ground_cond*T_ground);
    term1_w = W_zone_gain + air_density*(mixing*m_v_zone + infiltration.*m_v_air) +  m_sys.*m_v_supply;%flows in m^3/s converted to kg/s, everything into kg/s of water 
    term1 = [term1_z;term1_s;term1_w];
    new = ([T_zone;T_surf;m_v_zone] - term1./den).*exp(-den./Cap*dt/sub_steps) + term1./den;
    T_zone = new(1:z);
    T_surf = new(z+1:z+s);
    m_v_zone = new(z+s+1:end);
end

% Q_surfaces = (h_A_T - h_A.*T_new');
% Q_mixing = air_density*spec_heat_zone.*(mixing*T_setpoint - sum(mixing,2).*T_new');%heat into zone from another zone
% Q_infiltration = air_density*infiltration.*spec_heat_air.*(Tair_z - T_new');
% Q_hvac = m_sys.*spec_heat_supply.*(T_supply-T_new');

% T_interior = building.ctf.map_sur2zone'*T_zone;
% T_surf_interior = T_new(building.ctf.sur_state1)';
% T_surf_exterior = T_new(building.ctf.sur_state2)';
% Q_int_cond = (building.ctf.R_pos + building.ctf.R_neg)*T_surf; %sum should be zero
% Q_from_zone = (Q_surf_absorb(building.ctf.sur_state1) + h.*(T_interior - T_surf_interior).*building.surfaces.area + [T_surf_interior;T_windows]'*h_rad - sum(h_rad,1)'.*T_surf_interior); %at steady state total should equal - total(Q_from_environment)
% Q_from_environment = (Q_surf_absorb(building.ctf.sur_state2) + (h_s.*(T_exterior-T_surf_exterior) + h_sky.*(T_sky-T_surf_exterior)).*building.ctf.exterior_area + ground_cond.*(T_ground - T_surf_exterior));
% Q_in_sur =  Q_int_cond + building.ctf.interior_surface*Q_from_zone + building.ctf.exterior_surface*Q_from_environment;
% net_2_surf = sum(Q_in_sur );
% net_rad_imbalance = sum(h_rad*T_surf_interior - sum(h_rad,2).*T_surf_interior);
end