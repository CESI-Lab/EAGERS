function [loads,gains,occupancy,mixing,infiltration,T_set,w_set,vect_2_sun,frost,T_mains,water_heat] = ...
    zone_loads(building,date,profile,weather,T_zone,m_v_zone,dt,frost)
m_v_air = psychometric('P',weather.PressurePa,'Tdp',weather.DewpointC,'w');
T_air = weather.DrybulbC*ones(1,length(building.ctf.z_height)) - ones(length(date),1)*0.0065*building.ctf.z_height';%outdoor air temp next to each zone
n_s = length(date);
z = length(building.zones.name);
if isempty(T_zone)
    T_zone = 21*ones(n_s,z);%Assumed zone air temperature (ideally a function as building operates)
end
T_mains = water_mains(building.water_main_temperature,date);
[T_set,w_set] = thermostat_humidistat(building,profile,T_zone);%% HVAC temperature Schedules

%% Occupancy and latent heat from occupants
cat = {'convected';'radiant';'latent';'visible';'lost'};
for j = 1:1:length(cat)
    loads.(cat{j}) = zeros(n_s,z);
end

occupancy = zeros(n_s,z);
loads.occupancy = zeros(n_s,z);
loads.internal_gain = zeros(n_s,z);
n = length(building.occupancy.zone);
for i = 1:1:n
    z_i = building.occupancy.zone(i);
    sched_i = building.occupancy.schedule{i};
    new_occ = profile.(sched_i)*building.occupancy.nominal(i)*building.zones.floor_area(z_i);
    occupancy(:,z_i) = occupancy(:,z_i) + new_occ;
    if ~isempty(building.occupancy.activity{i})
        metabolic_rate = profile.(building.occupancy.activity{i});% W/person
    else
        metabolic_rate = 120;%120W/person
    end
    sensible_heat = 6.461927 + 0.946892*metabolic_rate + .0000255737*metabolic_rate.^2 + 7.139322*T_zone(:,z_i) - ...
        0.0627909*T_zone(:,z_i).*metabolic_rate + 0.0000589172*T_zone(:,z_i).*metabolic_rate.^2 - 0.198550*T_zone(:,z_i).^2 + ...
        0.000940018*T_zone(:,z_i).^2.*metabolic_rate - 0.00000149532*T_zone(:,z_i).^2.*metabolic_rate.^2;
    %%What do i do with clothing & work efficiency?
    loads.radiant(:,z_i) = loads.radiant(:,z_i) + building.occupancy.radiant(i)*new_occ.*sensible_heat;%radiant is a portion of the sensible heat
    loads.convected(:,z_i) = loads.convected(:,z_i) + building.occupancy.convected(i)*new_occ.*sensible_heat;
    loads.latent(:,z_i) = loads.latent(:,z_i) + new_occ.*(metabolic_rate - sensible_heat);
    loads.occupancy(:,z_i) = loads.occupancy(:,z_i) + new_occ.*metabolic_rate;
end
loads.internal_gain = loads.internal_gain + loads.occupancy;

loads.multiplier = building.zones.multiplier;
loads.area = building.zones.floor_area;

%% Infiltration and Mixing
infiltration = zeros(n_s,z);%volumetric flows (m3/s)
for j = 1:1:length(building.infiltration.nominal)
    infiltration(:,building.infiltration.zone(j)) = building.infiltration.nominal(j)*profile.(building.infiltration.schedule{j});%volumetric flows (m3/s)
end
mixing = zeros(z,z,n_s);%volumetric flows (m3/s)
for j = 1:1:length(building.mixing.nominal)
    if strcmpi(building.mixing.type{j},'Flow/Person')
        mixing(building.mixing.receiving_zone(j),building.mixing.source_zone(j),:) = building.mixing.nominal(j)*occupancy(:,building.mixing.source_zone(j)).*profile.(building.mixing.schedule{j});%volumetric flows (m3/s)
    else
        mixing(building.mixing.receiving_zone(j),building.mixing.source_zone(j),:) = building.mixing.nominal(j)*profile.(building.mixing.schedule{j});%volumetric flows (m3/s)
    end
end

%% Equipment and Lighting
prop = {'lighting_internal';'plug_load';'gas_load';};
for j = 1:1:length(prop)
    loads.(prop{j}) = zeros(n_s,z);
    if ~isempty(building.(prop{j}))
        n = length(building.(prop{j}).zone);
        for i = 1:1:n
            z_i = building.(prop{j}).zone(i);
            sched_i = building.(prop{j}).schedule{i};
            loads.(prop{j})(:,z_i) = loads.(prop{j})(:,z_i) + profile.(sched_i)*building.(prop{j}).nominal(i)*building.zones.floor_area(z_i); 
            loads.internal_gain(:,z_i) = loads.internal_gain(:,z_i) + loads.(prop{j})(:,z_i);
            for c = 1:1:length(cat)
                if isfield(building.(prop{j}),cat{c})
                    loads.(cat{c})(:,z_i) = loads.(cat{c})(:,z_i) + loads.(prop{j})(:,z_i)*building.(prop{j}).(cat{c})(i);
                end
            end
        end
    end
end

%% Solar
[sunrise, sunset, azimuth, zenith] = solar_calc(building.location.Longitude,building.location.Latitude,building.location.TimeZone,date);
vect_2_sun = [sind(zenith).*sind(azimuth), sind(zenith).*cosd(azimuth), cosd(zenith)];%% normal vector pointed at sun from azimuth (compass heading of sun) and zenith (angle beween sun and normal to surface of earth)
ext_gain = zeros(length(date),length(building.surfaces.exterior.cos_phi));%solar gain per unit area of exterior
if any((weather.DNIWm2+weather.DHIWm2)>0)
    F_s = 0.5*(1+building.surfaces.exterior.cos_phi)';
    F_g = 0.5*(1-building.surfaces.exterior.cos_phi)';
    ext_gain = weather.DHIWm2*(building.surfaces.exterior.solar_absorptance'.*F_s) + building.site.ground_reflect*weather.DHIWm2*(building.surfaces.exterior.solar_absorptance'.*F_g);
    cos_theta = vect_2_sun(:,1)*building.surfaces.exterior.normal(:,1)' + vect_2_sun(:,2)*building.surfaces.exterior.normal(:,2)' + vect_2_sun(:,3)*building.surfaces.exterior.normal(:,3)';%dot product of two vectors
    SF = 1; %shading factor, needs development, pg 211
    ext_gain = ext_gain +  weather.DNIWm2*building.surfaces.exterior.solar_absorptance'*SF.*cos_theta.*(cos_theta>0);%Absorbed direct normal radiation
    g = f_index('Ground',building.surfaces.exterior.boundary);
    ext_gain(:,g) = 0;
end
        
%% Refrigeration
[loads.case,ref_gain,frost] = refrigerated_cases(building,profile,T_zone,m_v_zone,frost,dt);
[loads.rack,rack_gain] = compressor_racks(building,loads.case.evap,T_air(:,1),T_zone,m_v_air);

%% Water Equipment
[water_gain, water_heat] = water_equipment_gains(building,profile,T_zone,m_v_zone,T_mains,dt);

%% Exterior loads
frac_dark = astronomical_clock_controlled_exterior_lights(date,sunrise,sunset);
loads.exterior.lighting = zeros(length(date),1);
loads.exterior.equipment = zeros(length(date),1);
for j = 1:1:length(building.exterior.name)
    prof = profile.(building.exterior.schedule{j});
    if strcmp(building.exterior.type{j},'Lights')
        loads.exterior.lighting = loads.exterior.lighting + prof*building.exterior.nominal(j).*frac_dark;
    else
        loads.exterior.equipment = loads.exterior.equipment + prof*building.exterior.nominal(j);
    end
end

%% grouping terms
loads.latent = loads.latent + water_gain.latent;
gains.zone_sensible = loads.convected + rack_gain.zone + ref_gain.zone_sensible + water_gain.sensible;
gains.zone_latent = loads.latent + ref_gain.zone_latent;
gains.interior_surface = (loads.radiant + loads.visible)*building.ctf.map_zone2sur;
gains.exterior_surface = ext_gain;
gains.hvac_sensible = rack_gain.hvac + ref_gain.hvac_sensible;
gains.hvac_latent = ref_gain.hvac_latent;
gains.window.visible = (loads.visible*building.ctf.map_zone2win)./(ones(n_s,1)*building.windows.area');
gains.window.long_wave = (loads.radiant*building.ctf.map_zone2win)./(ones(n_s,1)*building.windows.area');
gains.ref_rack_water_sensible = rack_gain.zone + ref_gain.zone_sensible + water_gain.sensible;
gains.ref_rack_water_latent = water_gain.latent + ref_gain.zone_latent;
end%Ends function zone_loads

function [T_set,w_set] = thermostat_humidistat(building,profile,T_zone)
%% HVAC temperature Schedules
T_set.heat = T_zone;
T_set.cool = T_zone;
w_set.humidify = zeros(length(T_zone),1);
w_set.dehumidify = ones(length(T_zone),1);
T_set.no_hvac = true(1,length(building.zones.name));
for i = 1:1:length(building.zone_controls.name)
    z = f_index(building.zone_controls.zone{i},building.zones.name);
    switch building.zone_controls.type{i}
        case 'Thermostat'
            T_set.no_hvac(z) = false;
            type = profile.(building.zone_controls.control_type_schedule{i});
            c_avail = building.zone_controls.control_type{i};
            c_type = zeros(length(c_avail),1);
            for j = 1:1:length(c_avail)
                switch c_avail{j}
                    case 'ThermostatSetpoint:SingleHeating'
                        c_type(j) = 1;
                    case 'ThermostatSetpoint:SingleCooling'
                        c_type(j) = 2;
                    case 'ThermostatSetpoint:SingleHeatingOrCooling'
                        c_type(j) = 3;
                    case 'ThermostatSetpoint:DualSetpoint'
                        c_type(j) = 4;
                end
            end
            k = f_index(type(1),c_type);%need to generalize to work when asking for a time series of profile
            t = f_index(building.zone_controls.control_name{i}{k},building.thermostat.name);
            switch c_avail{k}
                case 'ThermostatSetpoint:SingleHeating'
                    T_set.heat(:,z) = profile.(building.thermostat.heating{t});
                case 'ThermostatSetpoint:SingleCooling'
                    T_set.cool(:,z) = profile.(building.thermostat.cooling{t});
                case 'ThermostatSetpoint:SingleHeatingOrCooling'
                    T_set.heat(:,z) = profile.(building.thermostat.heating{t});
                    T_set.cool(:,z) = profile.(building.thermostat.heating{t});
                case 'ThermostatSetpoint:DualSetpoint'
                    T_set.heat(:,z) = profile.(building.thermostat.heating{t});
                    T_set.cool(:,z) = profile.(building.thermostat.cooling{t});
            end
        case 'Humidistat'
            w_set.humidify(:,z) = profile.(building.zone_controls.humidify_schedule{i}{1});
            w_set.dehumidify(:,z) = profile.(building.zone_controls.dehumidify_schedule{i}{1});
    end
end
end%Ends function thermostat

function [loads,ref_gain,frost] = refrigerated_cases(building,profile,T_zone,m_v_zone,frost,dt)
%%See page 1343 of the reference manual
bc = building.cases;
n_s = length(T_zone(:,1));
n_c = length(bc.name);
loads.electric = zeros(n_s,n_c);
loads.evap = zeros(n_s,n_c);
ref_gain.zone_sensible = 0*T_zone; %cumulative sensible heat to all refrigerated cases from the building zone
ref_gain.zone_latent = 0*T_zone;
ref_gain.hvac_sensible = 0*T_zone; %cumulative sensible heat credits applied to HVAC
ref_gain.hvac_latent = 0*T_zone;
for i = 1:1:n_c
    z = f_index(bc.zone{i},building.zones.name);
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
            if length(dt) == 1
                [frost(i),Q_def] = frost_calc(frost(i),new_frost,P_def,max_melt_frost,dt);
            else
                for t = 1:1:length(dt)
                    if t == 1
                        [frost(t,i),Q_def(t)] = frost_calc(0,new_frost(t),P_def(t),max_melt_frost(t),dt(t));
                    else
                        [frost(t,i),Q_def(t)] = frost_calc(frost(t-1,i),new_frost(t),P_def(t),max_melt_frost(t),dt(t));
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
    RH = psychometric('Tdb',T_zone(z),'w',m_v_zone(z),'RH');
    Tdp = psychometric('Tdb',bc.temperature(i),'w',m_v_zone(z),'Tdp');

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
        Tdp_rated = psychometric('Tdb',bc.rated_ambient_T(i),'RH',bc.rated_ambient_RH(i),'Tdp');
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
    ref_gain_sensible = (P_lights.*(1-bc.light_to_case(i)) + P_as.*(1-bc.anti_sweat_heat_to_case(i)) - Q_sensible_credit).*profile.(sched_i);
    loads.evap(:,i) = (Q_sensible_credit + P_fan + P_lights*bc.light_to_case(i) + P_as*bc.anti_sweat_heat_to_case(i) + Q_restock + Q_def).*profile.(sched_i);
    %% Latent loads
    if strcmpi(bc.latent_credit_curve_type{i},'CaseTemperatureMethod')
        latent_ratio = 1 - (bc.rated_ambient_RH(i) - RH)*eval_curve(building.curves,bc.latent_credit_curve_name{i},bc.temperature(i));
    elseif strcmpi(bc.latent_credit_curve_type{i},'RelativeHumidityMethod')
        latent_ratio = eval_curve(building.curves,bc.latent_credit_curve_name{i},RH);
    elseif strcmpi(bc.latent_credit_curve_type{i},'DewpointMethod')
        latent_ratio = eval_curve(building.curves,bc.latent_credit_curve_name{i},Tdp);
    else
        disp('check names of case latent credit curve')
    end
    if latent_ratio>1 
        latent_ratio = 1;
    elseif latent_ratio<0
        latent_ratio = 0;
    end
    ref_gain_latent = bc.capacity_per_length(i)*L*bc.latent_heat_ratio(i)*bc.runtime_fraction(i)*credit_frac*latent_ratio;

    %%split into zone & hvac due to return_air fraction
    ref_gain.hvac_sensible(:,z) = ref_gain.hvac_sensible(:,z) + ref_gain_sensible*RAF;
    ref_gain.zone_sensible(:,z) = ref_gain.zone_sensible(:,z) + ref_gain_sensible*(1-RAF);
    ref_gain.hvac_latent(:,z) = ref_gain.hvac_latent(:,z) + ref_gain_latent*RAF;
    ref_gain.zone_latent(:,z) = ref_gain.zone_latent(:,z) + ref_gain_latent*(1-RAF);

    loads.electric(:,i) = (P_fan + P_lights + P_as).*profile.(sched_i);
    if strcmpi(bc.defrost_type{i},'electric')
        loads.electric(:,i) = loads.electric(:,i) + P_def.*profile.(sched_i);
    end
end
end%Ends function refrigerated_cases

function [loads,rack_gain] = compressor_racks(building,case_evap_loads,T_air,T_zone,m_v_air)
n_s = length(T_air);
n_r = length(building.racks.name);
loads.electric = zeros(n_s,n_r);
loads.water = zeros(n_s,n_r);
loads.evap = zeros(n_s,n_r);
Q_reclaim = zeros(n_s,n_r);
rack_gain.zone = zeros(n_s,length(building.zones.name));
rack_gain.hvac = zeros(n_s,length(building.zones.name));
for i = 1:1:n_r
    loads.evap(:,i) = sum(case_evap_loads(:,building.racks.cases{i}),2);
    if strcmpi(building.racks.heat_reject_location{i},'outdoors')
        z = 0;
        if strcmpi(building.racks.condensor_type{i},'aircooled')
            T_eval = T_air;
        elseif strcmpi(building.racks.condensor_type{i},'evapcooled')
            Twb = psychometric('Tdb',T_air,'w',m_v_air,'Twb');
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
    COP = building.racks.design_COP(i)*eval_curve(building.curves,building.racks.COP_curve{i},T_eval);
    if strcmpi(building.racks.condensor_type{i},'watercooled')
        fan_power = 0;
        loads.water(:,i) = loads.evap(:,i)./COP;
    elseif isempty(building.racks.fan_power_temperature_curve{i})
        fan_power = building.racks.fan_power(i);
    else
        fan_power = building.racks.fan_power(i)*eval_curve(building.curves,building.racks.fan_power_temperature_curve{i},T_eval);
    end
    if ~strcmpi(building.racks.heat_reject_location{i},'outdoors') || strcmpi(building.racks.condensor_type{i},'aircooled') || strcmpi(building.racks.condensor_type{i},'evapcooled')
        loads.electric(:,i) = loads.evap(:,i)./COP + fan_power;
    end
    if z>0
        %if any are walk-in
%         Q_rack_2_zone(:,z) = loads.evap(:,i)/COP + fan_power;
%         else
        rack_gain.zone(:,z) = sum(case_evap_loads(:,building.racks.cases{i}).*(1-building.cases.return_air_frac(building.racks.cases{i})),2)/sum(case_evap_loads(:,building.racks.cases{i}),2)*(loads.evap(:,i)/COP + fan_power);
        rack_gain.hvac(:,z) = (loads.evap(:,i)./COP + fan_power) - rack_gain.zone(:,z);
    else
        %possible to reclaim up to 30% waste heat
        Q_reclaim(:,i) = 0.3*loads.evap(:,i).*(1 + 1./COP);
    end
end
end%Ends function compressor_racks

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

function [frost,Q_def] = frost_calc(frost,new_frost,P_def,max_melt_frost,dt)
frost = frost + new_frost;
Q_def = 0;
if P_def>0
    if frost>max_melt_frost
        frost = frost-max_melt_frost;
    else
        frost = 0;
        Q_def = P_def - frost*335000./dt;
    end 
end
end%ends function frost_calc

function T_mains = water_mains(main,date)
if isempty(main)
    T_mains = 15;
else
    [year,~,~] = datevec(date);
    julian_day = ceil(date - datenum(year,1,1));
    T_out_avg = main.outdoor_average*9/5+32; %convert to F 
    T_max_diff = main.max_temp_difference*9/5;%convert to delta F
    ratio = 0.4 + 0.01*(T_out_avg-44);
    lag = 35 - 1.0*(T_out_avg-44);
    T_mains = (T_out_avg+6) + ratio*(T_max_diff/2)*sin((0.986*(julian_day - 15 - lag) -90)*pi/180);
    T_mains = (T_mains-32)*5/9;%convert back to C
end
end%Ends function water_mains

function [water_gain, water_heat] = water_equipment_gains(building,profile,T_zone,m_v_zone,T_mains,dt)
Cp_water = 4186; %J/kg*K
den_water = 995;%kg/m^3
air_density = 1.204; %kg/m^3 %standard air density (1.204 kg/m3) adjusted for the local barometric pressure (standard barometric pressure corrected for altitude, ASHRAE 1997 HOF pg. 6.1).
water_gain.sensible = 0*T_zone;
water_gain.latent = 0*T_zone;
water_heat = zeros(length(dt),1);
for j = 1:1:length(building.water_use.name)
    if strcmp(building.water_use.type{j},'Equipment')
        [~,~,T_target,T_hot,T_cold] = water_equipment(building,profile,T_mains,[],j);
        %%find gain into zones
        z = f_index(building.water_use.zone(j),building.zones.name);
        m_mix = building.water_use.peak_flow(j)*profile.(building.water_use.flow_schedule{j})*den_water; 
        m_hot = m_mix.*(T_target- T_cold)./(T_hot - T_cold);
        Q_total = m_mix*Cp_water.*(T_target - T_zone(:,z));%Watts of energy relative to zone
%         T_wb = psychometric('Tdb',T_zone(:,z),'w',m_v_zone(:,z),'Twb');%wet-bulb temperature     
%         m_v_sat = psychometric('Twb',T_wb,'w');%mass fraction of water in saturated air     
        m_v_sat = psychometric('Twb',T_zone(:,z),'w');%mass fraction of water in saturated air     
        e_num = f_index(building.water_use.name{j},building.plant_demand_equip.name);
        loop_connect = building.plant_demand_equip.loop(e_num);
        if ~any(strcmpi(building.plant_loop.type(loop_connect),'Heating'))
            water_heat = water_heat + (m_hot*Cp_water*(T_hot - T_cold))*building.zones.multiplier(z);%only if not on a plant loop
        end
        water_gain.latent(:,z) = profile.(building.water_use.latent_frac_schedule{j}).*Q_total;
%         water_gain.latent(:,z) = min(profile.(building.water_use.latent_frac_schedule{j}).*Q_total,(m_v_sat - m_v_zone(:,z))*air_density*building.zones.volume(z)*2260000./dt);
        water_gain.sensible(:,z) = Q_total.*profile.(building.water_use.sensible_frac_schedule{j});%energy added to zone (without zone multiplier)
    end
end
end%Ends function water_equipment_gains