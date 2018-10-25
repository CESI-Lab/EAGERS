function [net_building,loads,occupancy,mixing,infiltration,Q_zone_gain,Q_surf_absorb,...
    W_zone_gain,internal_irrad_window,T_set,cos_phi_windows,frost] = ...
    zone_loads(building,date,profile,solar_normal,solar_diffuse,T_zone,m_v_zone,dt,T_air,m_v_air,frost)

n_s = length(date);
z = length(building.zones.name);
if isempty(T_zone)
    T_zone = 21*ones(n_s,z);%Assumed zone air temperature (ideally a function as building operates)
end
% Schedules
sched = fieldnames(building.schedule);
T_set.heat = T_zone;
T_set.cool = T_zone;
T_set.supply_h = nan(1,z);
T_set.supply_c = nan(1,z);
T_set.no_hvac = false(1,z);
infiltration = zeros(n_s,z);
mixing = zeros(z,z,n_s);
for i = 1:1:z
    if ~isempty(building.zones.thermostat{i,1})
        T_set.heat(:,i) = profile.(building.zones.thermostat{i,1});
        T_set.cool(:,i) = profile.(building.zones.thermostat{i,2});
        T_set.supply_h(i) = building.zones.t_supply_h(i);
        T_set.supply_c(i) = building.zones.t_supply_c(i);
    else
        T_set.no_hvac(i) = true;
    end
    if building.zones.infiltration.nominal(i)>0
        infiltration(:,i) = building.zones.infiltration.nominal(i)*profile.(building.zones.infiltration.schedule{i});
    end
    for j = 1:1:length(building.zones.mixing.source{i})
        k = building.zones.mixing.source{i}(j);
        mixing(i,k,:) = building.zones.mixing.nominal{i}(j)*profile.(building.zones.mixing.schedule{i}{j});
    end
end

cat = {'convected';'radiant';'latent';'visible';'lost'};
for j = 1:1:length(cat)
    loads.(cat{j}) = zeros(n_s,z);
end

occupancy = zeros(n_s,z);
loads.occupancy = zeros(n_s,z);
loads.internal_gain = zeros(n_s,z);
n = length(building.gains.occupancy.zone);
for i = 1:1:n
    z_i = building.gains.occupancy.zone(i);
    sched_i = building.gains.occupancy.schedule{i};
    new_occ = profile.(sched_i)*building.gains.occupancy.nominal(i)*building.zones.floor_area(z_i);
    occupancy(:,z_i) = occupancy(:,z_i) + new_occ;
    if ~isempty(building.gains.occupancy.activity{i})
        metabolic_rate = profile.(building.gains.occupancy.activity{i});% W/person
    else
        metabolic_rate = 120;%120W/person
    end
    sensible_heat = 6.461927 + 0.946892*metabolic_rate + .0000255737*metabolic_rate.^2 + 7.139322*T_zone(:,z_i) - ...
        0.0627909*T_zone(:,z_i).*metabolic_rate + 0.0000589172*T_zone(:,z_i).*metabolic_rate.^2 - 0.198550*T_zone(:,z_i).^2 + ...
        0.000940018*T_zone(:,z_i).^2.*metabolic_rate - 0.00000149532*T_zone(:,z_i).^2.*metabolic_rate.^2;
    %%What do i do with clothing & work efficiency?
    loads.radiant(:,z_i) = loads.radiant(:,z_i) + building.gains.occupancy.radiant(i)*new_occ.*sensible_heat;%radiant is a portion of the sensible heat
    loads.convected(:,z_i) = loads.convected(:,z_i) + building.gains.occupancy.convected(i)*new_occ.*sensible_heat;
    loads.latent(:,z_i) = loads.latent(:,z_i) + new_occ.*(metabolic_rate - sensible_heat);
    loads.occupancy(:,z_i) = loads.occupancy(:,z_i) + new_occ.*metabolic_rate;
end
loads.internal_gain = loads.internal_gain + loads.occupancy;

loads.multiplier = building.zones.multiplier;
loads.area = building.zones.floor_area;


prop = {'lighting_internal';'plug_load';'gas_load';};
for j = 1:1:length(prop)
    loads.(prop{j}) = zeros(n_s,z);
    if ~isempty(building.gains.(prop{j}))
        n = length(building.gains.(prop{j}).zone);
        for i = 1:1:n
            z_i = building.gains.(prop{j}).zone(i);
            sched_i = building.gains.(prop{j}).schedule{i};
            loads.(prop{j})(:,z_i) = loads.(prop{j})(:,z_i) + profile.(sched_i)*building.gains.(prop{j}).nominal(i)*building.zones.floor_area(z_i); 
            loads.internal_gain(:,z_i) = loads.internal_gain(:,z_i) + loads.(prop{j})(:,z_i);
            for c = 1:1:length(cat)
                if isfield(building.gains.(prop{j}),cat{c})
                    loads.(cat{c})(:,z_i) = loads.(cat{c})(:,z_i) + loads.(prop{j})(:,z_i)*building.gains.(prop{j}).(cat{c})(i);
                end
            end
        end
    end
end

[sunrise, sunset, azimuth, zenith] = solar_calc(building.location.Longitude,building.location.Latitude,building.location.TimeZone,date);
vect_2_sun = [sind(zenith).*sind(azimuth), sind(zenith).*cosd(azimuth), cosd(zenith)];%% normal vector pointed at sun from azimuth (compass heading of sun) and zenith (angle beween sun and normal to surface of earth)
cos_phi_sun = max(0,vect_2_sun*building.surfaces.normal');%magnitde of both vectors is one so a dot b/|a||b| is just a dot b, %% If dot product of surface normal and vector pointed at sun is <0, set cos_phi = 0 because it is shading itself.
cos_phi_windows = cos_phi_sun(:,building.windows.surface);

[net_case_electric_loads,case_evap_loads,Q_ref_2_zone,Q_ref_2_HVAC,frost] = refrigerated_cases(building.cases,profile,building.HVAC.curves,T_zone,m_v_zone,frost,dt);
[net_rack_electric_loads,Q_rack_2_zone,Q_rack_2_HVAC] = compressor_racks(building,case_evap_loads,T_air(:,1),T_zone,m_v_air);
[hot_water_use,cold_water_use,Q_water_2_zone] = water_equipment(building.water_use,building.zones.volume,profile,T_zone,m_v_zone,dt);

Q_zone_gain = loads.convected + Q_rack_2_zone + Q_ref_2_zone.sensible + Q_water_2_zone.sensible;
Q_surf_absorb_interior = (loads.radiant + loads.visible)*building.ctf.map_zone2sur;
Q_surf_absorb_exterior = surf_solar_gain(building,solar_normal,solar_diffuse,cos_phi_sun);%solar gains
Q_surf_absorb = Q_surf_absorb_interior*building.ctf.interior_surface' + Q_surf_absorb_exterior*building.ctf.exterior_surface';

W_zone_gain = (loads.latent + Q_water_2_zone.latent + Q_ref_2_zone.latent)/2264705;%convert J/s to kg/s of water
loads.internal_radiation.visible = (loads.visible*building.ctf.map_zone2win)./(ones(n_s,1)*building.windows.area');
loads.internal_radiation.long_wave = (loads.radiant*building.ctf.map_zone2win)./(ones(n_s,1)*building.windows.area');
internal_irrad_window = loads.internal_radiation;

Q_HVAC_gain = Q_rack_2_HVAC + Q_ref_2_HVAC.sensible;
W_HVAC_gain = Q_ref_2_HVAC.latent;

frac_dark = astronomical_clock_controlled_exterior_lights(date,sunrise,sunset);
loads.exterior.lighting = zeros(length(date),1);
loads.exterior.equipment = zeros(length(date),1);
for j = 1:1:length(building.exterior.name)
    sched_i = sched{nonzeros((1:length(sched))'.*strcmpi(building.exterior.schedule{j},sched))};
    if strcmp(building.exterior.type{j},'Lights')
        loads.exterior.lighting = loads.exterior.lighting + profile.(sched_i)*building.exterior.nominal(j).*frac_dark;
    else
        loads.exterior.equipment = loads.exterior.equipment + profile.(sched_i)*building.exterior.nominal(j);
    end
end

net_building.interior_lighting = loads.lighting_internal*loads.multiplier;
net_building.equipment = loads.plug_load*loads.multiplier;
net_building.internal_gain = loads.internal_gain*loads.multiplier;
net_building.non_hvac_electric = net_building.interior_lighting + net_building.equipment + loads.exterior.lighting + loads.exterior.equipment + net_case_electric_loads + net_rack_electric_loads;
end%Ends function zone_loads

function [net_case_electric_loads,case_evap_loads,Q_ref_2_zone,Q_ref_2_HVAC,frost] = refrigerated_cases(bc,profile,curves,T_zone,m_v_zone,frost,dt)
%%See page 1343 of the reference manual
n_s = length(T_zone(:,1));
net_case_electric_loads = zeros(n_s,1);
Q_ref_2_zone.sensible = 0*T_zone; %cumulative sensible heat to all refrigerated cases from the building zone
Q_ref_2_zone.latent = 0*T_zone;
Q_ref_2_HVAC.sensible = 0*T_zone; %cumulative sensible heat credits applied to HVAC
Q_ref_2_HVAC.latent = 0*T_zone;
case_evap_loads = 0;
if ~isempty(bc)
    c = length(bc.name);
    case_evap_loads = zeros(n_s,c);
    for i = 1:1:c
        z = bc.zone(i);
        L = bc.length(i);
        RAF = bc.return_air_frac(i);
        sched_i = bc.schedule{i};
        Q_def = zeros(n_s,1);
        if strcmpi(bc.defrost_type{i},'none')
                P_fan = bc.fan_power_per_length(i)*L;
        elseif strcmpi(bc.defrost_type{i},'electric')
                sched_def = profile.(bc.defrost_schedule{i});%fraction of time defrosting
                sched_drip = profile.(bc.defrost_drip_down_schedule{i});
                P_def = bc.defrost_power_per_length(i)*L*sched_def;
                latent_ratio = 1;
                new_frost = bc.capacity_per_length(i)*L*bc.runtime_fraction(i)*bc.latent_heat_ratio(i)*latent_ratio*dt.*(1-sched_drip)/(335000 + 2498000);
                max_melt_frost = P_def.*dt/335000;
                for t = 1:1:length(dt)
                    if t == 1
                        frost(t,i) = new_frost(t);
                    else
                        frost(t,i) = frost(t-1,i) + new_frost(t);
                    end
                    if P_def(t)>0
                        if frost(t,i)>max_melt_frost(t)
                            frost(t,i) =frost(t,i)-max_melt_frost(t);
                        else
                            frost(t,i) = 0;
                            Q_def(t) = P_def(t) - frost(t,i)*335000./dt(t);
                        end 
                    end
                end
                P_fan = bc.fan_power_per_length(i)*L*(1-sched_def);
        elseif strcmpi(bc.defrost_type{i},'offcycle')
%                 sched_drip = profile.(bc.defrost_drip_down_schedule{i});
                P_fan = bc.fan_power_per_length(i)*L;
        elseif strcmpi(bc.defrost_type{i},'hotgas')
            
        elseif strcmpi(bc.defrost_type{i},'hotbrine')

        end
        P_H2O = 101.325*m_v_zone(z)./(0.62198 + m_v_zone(z));
        Tdb_K = T_zone(z)+273.15; %Tdb (Kelvin)
        satP = exp((-5.8002206e3)./Tdb_K + 1.3914993 - 4.8640239e-2*Tdb_K + 4.1764768e-5*Tdb_K.^2 - 1.4452093e-8*Tdb_K.^3 + 6.5459673*log(Tdb_K))/1000; %saturated water vapor pressure ASHRAE 2013 fundamentals eq. 6 in kPa valid for 0 to 200C
        RH = P_H2O./satP*100;%Relative humidity in percent
        Tdp = psychometric_Tdp(bc.temperature(i),m_v_zone(z));
            
        P_lights = bc.standard_lighting_per_unit_length(i)*L*profile.(bc.light_schedule{i});
        if strcmpi(bc.anti_sweat_control{i},'none')
            P_as = 0;
        elseif strcmpi(bc.anti_sweat_control{i},'constant')
            P_as = bc.anti_sweat_heater_per_length(i)*L;
        elseif strcmpi(bc.anti_sweat_control{i},'RelativeHumid') || strcmpi(bc.anti_sweat_control{i},'linear')
            P_as = bc.anti_sweat_heater_per_length(i)*L.*(1 - (bc.rated_ambient_RH(i) - RH)/(bc.rated_ambient_RH(i) - bc.humidity_at_zero_percent(i)));
        elseif strcmpi(bc.anti_sweat_control{i},'dewpoint')
            %%%
        elseif strcmpi(bc.anti_sweat_control{i},'heatbalance') || strcmpi(bc.anti_sweat_control{i},'HeatBalanceMethod')
            %% first find R_case at nominal conditions
            P = 101.325; % atmospheric pressure (kPa)
            Tdb_K = bc.rated_ambient_T(i)+273.15; %Tdb (Kelvin)
            satP = exp((-5.8002206e3)./Tdb_K + 1.3914993 - 4.8640239e-2*Tdb_K + 4.1764768e-5*Tdb_K.^2 - 1.4452093e-8*Tdb_K.^3 + 6.5459673*log(Tdb_K))/1000; %saturated water vapor pressure ASHRAE 2013 fundamentals eq. 6 in kPa valid for 0 to 200C
            P_H2O = bc.rated_ambient_RH(i)/100.*satP; % kPa
            m_v_rated = .621945*(P_H2O/(P-P_H2O));
            Tdp_rated = psychometric_Tdp(bc.rated_ambient_T(i),m_v_rated);
            R_case = (Tdp_rated - bc.temperature(i))./(bc.anti_sweat_heater_per_length(i)./bc.height(i) - (Tdp_rated - bc.rated_ambient_T(i))./0.3169);
            P_as = ((Tdp - T_zone(z)).*bc.height(i)/0.3169 + (Tdp - bc.temperature(i))./R_case).*L;
        end
        P_as = max(P_as,bc.minimum_anti_sweat_per_length(i)*L);
        Q_restock = profile.(bc.restock_schedule{i})*L;
        if length(bc.case_credit_fraction_schedule)>=i && ~isempty(bc.case_credit_fraction_schedule{i})
            credit_frac = profile.(bc.case_credit_fraction_schedule{i});
        else
            credit_frac = 1;
        end
        Q_sensible_credit_rated = (bc.capacity_per_length(i)*L*bc.runtime_fraction(i)*(1-bc.latent_heat_ratio(i)) - P_lights*bc.light_to_case(i) - P_as*bc.anti_sweat_heat_to_case(i) - P_fan);
        Q_sensible_credit = Q_sensible_credit_rated.*((T_zone(:,z) - bc.temperature(i))/(bc.rated_ambient_T(i) - bc.temperature(i))).*credit_frac;
        Q_ref_2_zone.sensible(:,z) = Q_ref_2_zone.sensible(:,z) + (P_lights.*(1-bc.light_to_case(i)) + P_as.*(1-bc.anti_sweat_heat_to_case(i)) - Q_sensible_credit).*profile.(sched_i);
        case_evap_loads(:,i) = (Q_sensible_credit + P_fan + P_lights*bc.light_to_case(i) + P_as*bc.anti_sweat_heat_to_case(i) + Q_restock + Q_def).*profile.(sched_i);
        %% Latent loads
        if strcmpi(bc.latent_credit_curve_type{i},'CaseTemperatureMethod')
            latent_ratio = 1 - (bc.rated_ambient_RH(i) - RH)*eval_curve(curves,bc.latent_credit_curve_name{i},bc.temperature(i));
        elseif strcmpi(bc.latent_credit_curve_type{i},'RelativeHumidityMethod')
            latent_ratio = eval_curve(curves,bc.latent_credit_curve_name{i},RH);
        elseif strcmpi(bc.latent_credit_curve_type{i},'DewpointMethod')
            latent_ratio = eval_curve(curves,bc.latent_credit_curve_name{i},Tdp);
        else
            disp('check names of case latent credit curve')
        end
        if latent_ratio>1 
            latent_ratio = 1;
        elseif latent_ratio<0
            latent_ratio = 0;
        end
        Q_ref_2_zone.latent(:,z) = bc.capacity_per_length(i)*L*bc.latent_heat_ratio(i)*bc.runtime_fraction(i)*credit_frac*latent_ratio;
        
        %%split into zone & hvac due to return_air fraction
        Q_ref_2_HVAC.sensible(:,z) = Q_ref_2_zone.sensible(:,z)*RAF;
        Q_ref_2_zone.sensible(:,z) = Q_ref_2_zone.sensible(:,z)*(1-RAF);
        Q_ref_2_HVAC.latent(:,z) = Q_ref_2_zone.latent(:,z)*RAF;
        Q_ref_2_zone.latent(:,z) = Q_ref_2_zone.latent(:,z)*(1-RAF);
        
        net_case_electric_loads = net_case_electric_loads + (P_fan + P_lights + P_as).*profile.(sched_i);
        if strcmpi(bc.defrost_type{i},'electric')
            net_case_electric_loads = net_case_electric_loads + P_def.*profile.(sched_i);
        end
    end
end
end%Ends function refrigerated_cases

function [net_rack_electric_loads,Q_rack_2_zone,Q_rack_2_HVAC] = compressor_racks(building,case_evap_loads,T_air,T_zone,m_v_air)
n_s = length(T_air);
n_r = length(building.racks.name);
rack_electric_loads = zeros(n_s,n_r);
rack_water_loads = zeros(n_s,n_r);
rack_evap_loads = zeros(n_s,n_r);
Q_reclaim = zeros(n_s,n_r);
Q_rack_2_zone = zeros(n_s,length(building.zones.name));
Q_rack_2_HVAC = zeros(n_s,length(building.zones.name));
for i = 1:1:n_r
    rack_evap_loads(:,i) = sum(case_evap_loads(:,building.racks.cases{i}),2);
    if strcmpi(building.racks.heat_reject_location{i},'outdoors')
        z = 0;
        if strcmpi(building.racks.condensor_type{i},'aircooled')
            T_eval = T_air;
        elseif strcmpi(building.racks.condensor_type{i},'evapcooled')
            Twb = psychometric_Twb(T_air,m_v_air);
            T_effective = T_wb + (1-evap_condenser_effectiveness(i))*(T_air - Twb);
            T_eval = T_effective;
        else % if strcmpi(building.racks.condensor_type{i},'watercooled') 
            %% need more details of water loop heat rejection
            T_eval = T_return_water;
        end
    else
        z = nonzeroes((1:length(building.zones.name))'.*strcmpi(building.racks.heat_reject_location{i},building.zones.name));
        T_eval = min(T_zone(z),[],2);
    end
    COP = building.racks.design_COP(i)*eval_curve(building.HVAC.curves,building.racks.COP_curve{i},T_eval);
    if strcmpi(building.racks.condensor_type{i},'watercooled')
        fan_power = 0;
        rack_water_loads(:,i) = rack_evap_loads(:,i)./COP;
    elseif isempty(building.racks.fan_power_temperature_curve{i})
        fan_power = building.racks.fan_power(i);
    else
        fan_power = building.racks.fan_power(i)*eval_curve(building.HVAC.curves,building.racks.fan_power_temperature_curve{i},T_eval);
    end
    if ~strcmpi(building.racks.heat_reject_location{i},'outdoors') || strcmpi(building.racks.condensor_type{i},'aircooled') || strcmpi(building.racks.condensor_type{i},'evapcooled')
        rack_electric_loads(:,i) = rack_evap_loads(:,i)./COP + fan_power;
    end
    if z>0
        %if any are walk-in
%         Q_rack_2_zone(:,z) = rack_evap_loads(:,i)/COP + fan_power;
%         else
        Q_rack_2_zone(:,z) = sum(case_evap_loads(:,building.racks.cases{i}).*(1-building.cases.return_air_frac(building.racks.cases{i})),2)/sum(case_evap_loads(:,building.racks.cases{i}),2)*(rack_evap_loads(:,i)/COP + fan_power);
        Q_rack_2_HVAC(:,z) = (rack_evap_loads(:,i)./COP + fan_power) - Q_rack_2_zone(:,z);
    else
        %possible to reclaim up to 30% waste heat
        Q_reclaim(:,i) = 0.3*rack_evap_loads(:,i).*(1 + 1./COP);
    end
end
net_rack_electric_loads = sum(rack_electric_loads,2);
end%Ends function compressor_racks

function [m_dot_drain,T_drain,Q_water_2_zone] = water_equipment(we,zone_volume,profile,T_zone,m_v_zone,dt)
Q_water_2_zone.sensible = 0*T_zone;
Q_water_2_zone.latent = 0*T_zone;
air_density = 1.225; %kg/m^3
if ~isempty(we)
    for i = 1:1:length(we.name)
        z = we.zone(i);
        T_sat = T_zone(:,z)+273;
        P = 101; % room pressure (kPa)
        P_H2O_sat = exp((-5.8002206e3)./T_sat + 1.3914993 - 4.8640239e-2*T_sat + 4.1764768e-5*T_sat.^2 - 1.4452093e-8*T_sat.^3 + 6.5459673*log(T_sat))/1000; %saturated water vapor pressure
        m_v_sat = .621945*(P_H2O_sat./(P-P_H2O_sat));%mass fraction of water in air 

        T_target = profile.(we.target_temperature_schedule{i});
        m_target = we.peak_flow(i)*profile.(we.flow_schedule{i})*995; % 995kg/m^3 at 30C
        m_dot_evap = profile.(we.latent_frac_schedule{i}).*min(m_target,(m_v_sat - m_v_zone(:,z))*air_density*zone_volume(z)./dt);
        Q_latent = m_dot_evap*2260000;
        Q_water_2_zone.latent(:,z) = Q_water_2_zone.latent(:,z) + Q_latent;
        Q_sensible = m_target*4186.*(T_target - T_zone(:,z)).*profile.(we.sensible_frac_schedule{i});
        Q_water_2_zone.sensible(:,z) = Q_water_2_zone.sensible(:,z) + Q_sensible;
        m_dot_drain = m_target - m_dot_evap;
        T_drain = (m_target.*4186.*T_target - Q_sensible - Q_latent)./(m_dot_drain*4186);
    end
else
    m_dot_drain = 0;
    T_drain = 0;
end
end%Ends function water_equipment

function surface_gain = surf_solar_gain(building,solar_normal,solar_diffuse,cos_phi)
%% Solar radiation hitting surfaces
n = length(building.surfaces.name);
exterior = strcmp(building.surfaces.boundary,'Outdoors');
cos_phi(:,~exterior) = 0;
solar_rad = cos_phi.*(solar_normal*ones(1,n)) + solar_diffuse*exterior';
%% gain onto surfaces
surfaces = building.surfaces;
n = length(surfaces.name);
surface_gain = zeros(length(solar_normal),length(building.ctf.exterior_area));
k = 0;
for i = 1:1:n
    if strcmp(surfaces.boundary{i},'Outdoors')%calculate energy absorbed
        mat = building.material.(building.construction.(surfaces.construct{i}){1});
        k = k+1;
        surface_gain(:,k) = building.ctf.exterior_area(k)*mat.solar_absorptance*solar_rad(:,i);
    end
end
end%Ends function surf_solar_gain

function frac_dark = astronomical_clock_controlled_exterior_lights(date,sunrise,sunset)
n_s = length(date);
shift = .04;%.125; % fraction of day exterior lights remain on after sunrise
frac_dark = ones(n_s,1);
if n_s == 1
    dt = [1,1];
else
    dt = 24*(date(2:end) - date(1:end-1));
    dt = [dt(1);dt];%shift things by 1 so frac_dark at first date can be calculated
end
date = [date(1)-dt(1);date];
hour =  24*(date - floor(date)); % hour of day
sun_up = hour(2:end)>(sunrise*24+shift) & hour(1:end-1)<=(sunrise*24+shift);%sunrise
s_up_frac = 1-((sunrise*24+shift) - hour(1:end-1))./dt;
frac_dark(sun_up) = s_up_frac(sun_up);
daytime = hour(2:end)<(sunset*24-shift);
frac_dark(daytime) = 0;
sun_down = hour(2:end)>(sunset*24-shift) & hour(1:end-1)<=(sunset*24-shift);%sunset
s_down_frac = ((sunset*24-shift) - hour(1:end-1))./dt;
frac_dark(sun_down) = s_down_frac(sun_down);
end%Ends function astronomical_clock_controlled_exterior_lights