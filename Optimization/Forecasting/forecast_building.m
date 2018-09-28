function [buildings,forecast] = forecast_building(forecast,date,buildings,options)
% This function estimates the energy profile of a building 
% date is a vector of points in datenum format at which you are requesting the electric cooling & heating power
% buildings is a structure of building parameters for each building
n_b = length(buildings);
n_s = length(date);

forecast.Building.ExternalGains = zeros(length(date),n_b);
for i = 1:1:n_b
    if isfield(buildings(i),'QPform')
        sgain = solar_gain(buildings(i),date,buildings(i).QPform.Location,forecast.Weather);
    else
        location.Longitude = 40;
        location.Latitude = 105;
        location.TimeZone = -6;
        sgain = solar_gain(buildings(i),date,location,forecast.Weather);
    end
    if strcmp(options.forecast,'Building')
        b_loads = building_loads(buildings(i),date,sgain);
        forecast.Building.InternalGains(:,i) = b_loads.InternalGains;
        forecast.Building.NonHVACelectric(:,i) = b_loads.Equipment + b_loads.InteriorLighting + b_loads.ExteriorLighting + b_loads.OtherLoads;
        if isfield(b_loads,'DCloads')
            forecast.Building.DCloads = b_loads.DCloads;
        end
    end
    forecast.Building.ExternalGains(:,i) = sgain.Walls + sgain.Roof;
end

if ~strcmp(options.solver,'NREL')
    forecast.Building.E0 = zeros(n_s,n_b);
    forecast.Building.C0 = zeros(n_s,n_b);
    forecast.Building.H0 = zeros(n_s,n_b);
    forecast.Building.Tmin = zeros(n_s,n_b);
    forecast.Building.Tmax = zeros(n_s,n_b);
    forecast.Building.Tset_H = zeros(n_s,n_b);
    forecast.Building.Tset_C = zeros(n_s,n_b);
    forecast.Building.Fan_Power = zeros(n_s,n_b);
    forecast.Building.AirFlow = zeros(n_s,n_b);
    forecast.Building.Tzone = zeros(n_s,n_b);
    forecast.Building.Twall = zeros(n_s,n_b);
    forecast.Building.Damper = zeros(n_s,n_b);
    for i = 1:1:n_b
        [cooling, heating, fan_power,zone,wall,damper] = building_profile(buildings(i),date,forecast.Building.InternalGains(:,i),forecast.Building.ExternalGains(:,i),forecast.Weather.DrybulbC,forecast.Weather.RHum,buildings(i).Tzone,buildings(i).Twall);
        forecast.Building.Tmin(:,i) = load_sched(buildings(i).Schedule,date,'TsetH') - 0.5*buildings(i).VariableStruct.Comfort;
        forecast.Building.Tmax(:,i) = load_sched(buildings(i).Schedule,date,'TsetC') + 0.5*buildings(i).VariableStruct.Comfort;
        [temperature_to_heat,temperature_to_cool] = equilib_temperature(buildings(i),date,forecast.Weather,zone(2:end),wall(2:end),cooling,heating,forecast.Building.InternalGains(:,i),forecast.Building.ExternalGains(:,i),forecast.Building.Tmin(:,i),forecast.Building.Tmax(:,i));
        forecast.Building.E0(:,i) = forecast.Building.NonHVACelectric(:,i) + fan_power(:,i);
        forecast.Building.C0(:,i) = cooling(:,i) ;
        forecast.Building.H0(:,i) = heating(:,i);
        if ~buildings(i).QPform.Heating
            %Can upgrade to include non-linear mapping of heating to electricity for built-in HVAC equipment
            forecast.Building.E0(:,i) = forecast.Building.E0(:,i) + (forecast.Building.H0(:,i)- buildings(i).QPform.UA*(zone(2:end)-temperature_to_heat))*buildings(i).QPform.H2E;
        end
        if ~buildings(i).QPform.Cooling
            %Can upgrade to include non-linear mapping of cooling to electricity for built-in HVAC equipment
            forecast.Building.E0(:,i) = forecast.Building.E0(:,i) + (forecast.Building.C0(:,i)- buildings(i).QPform.UA*(temperature_to_cool-zone(2:end)))*buildings(i).QPform.C2E;
        end
        forecast.Building.Tset_H(:,i) = temperature_to_heat;
        forecast.Building.Tset_C(:,i) = temperature_to_cool;
        forecast.Building.Fan_Power(:,i) = fan_power;
        forecast.Building.AirFlow(:,i) = fan_power/buildings(i).VariableStruct.FanPower;
        forecast.Building.Tzone(:,i) = zone(2:end);
        forecast.Building.Twall(:,i) = wall(2:end);
        forecast.Building.Damper(:,i) = damper;
    end
end
end%Ends function forecast_building

function [temperature_to_heat,temperature_to_cool] = equilib_temperature(build,date,weather,zone,wall,cooling,heating,internal_gains,external_gains,min_temp,max_temp)
rho_Air = 1.225; % kg/m^3
n = 20;
n_s = length(date);
temperature_to_heat = zeros(n_s,1);
temperature_to_cool = zeros(n_s,1);
if n_s == 1
    dt = 1;
else
    dt1 = date(2) - date(1);
    dt = (24*3600) * (date - [date(1)-dt1;date(1:end-1)]); % duration of each time segment [seconds]
end
%% Ambient dewpoint
pressure = 101.325; % atmospheric pressure (kPa)
dry_bulb_kelvin = weather.DrybulbC+273.15; %Tdb (Kelvin)
sat_press = exp((-5.8002206e3)./dry_bulb_kelvin + 1.3914993 - 4.8640239e-2*dry_bulb_kelvin + 4.1764768e-5*dry_bulb_kelvin.^2 - 1.4452093e-8*dry_bulb_kelvin.^3 + 6.5459673*log(dry_bulb_kelvin))/1000; %saturated water vapor pressure ASHRAE 2013 fundamentals eq. 6 in kPa valid for 0 to 200C
water_press = weather.RHum/100.*sat_press; % kPa
amb_spec_heat = 1.006 + 1.86*(.621945*(water_press./(pressure-water_press))); % kJ/kg*K
min_flow = build.VariableStruct.Volume*rho_Air*build.VariableStruct.AirChangePerHr*(1/3600); %kg/s of flow
for t = 1:1:n_s
    temperature_range = linspace(min_temp(t)-2*build.VariableStruct.Comfort,max_temp(t)+2*build.VariableStruct.Comfort,n)';
    net_gain = zone_thermal_load(build,wall(t),weather.DrybulbC(t),temperature_range,internal_gains(t),external_gains(t),dt(t));
    need_heat = (min(temperature_range,build.VariableStruct.ColdAirSet) - weather.DrybulbC(t))*amb_spec_heat(t)*min_flow*dt(t) - net_gain;
    if all(need_heat>0)
        temperature_to_heat(t) = temperature_range(1);
    elseif all(need_heat<0)
        temperature_to_heat(t) = temperature_range(end);
    else
        temperature_to_heat(t) = interp1(need_heat,temperature_range,0);
    end
    temperature_to_heat = max(temperature_to_heat,zone(t)-heating(t)*build.QPform.UA);
    need_cool = (weather.DrybulbC(t) - temperature_range)*amb_spec_heat(t)*min_flow*dt(t) + net_gain;
    if all(need_cool>0)
        temperature_to_cool(t) = temperature_range(end);
    elseif all(need_cool<0)
        temperature_to_cool(t) = temperature_range(1);
    else
        temperature_to_cool(t) = interp1(need_cool,temperature_range,0);
    end
    temperature_to_cool = min(temperature_to_cool,zone(t)+cooling(t)*build.QPform.UA);
end
end%Ends function equilib_temperature