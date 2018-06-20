function [buildings,b_forecast] = forecast_building(weather,date,buildings,prev_data,hist_prof,options)
% This function estimates the energy profile of a building 
% weather is an hourly weather profile (dry bulb, wet bulb, and relative humidity)
% date is a vector of points in datenum format at which you are requesting the electric cooling & heating power
% buildings is a structure of building parameters for each building
n_b = length(buildings);
n_s = length(date);
b_forecast.ExternalGains = zeros(length(date),n_b);
for i = 1:1:n_b
    sgain = solar_gain(buildings(i),date,buildings(i).QPform.Location,weather);
    b_forecast.ExternalGains(:,i) = sgain.Walls + sgain.Roof;
    if strcmp(options.forecast,'Building')
        b_loads = building_loads(buildings(i),date,sgain);
        b_forecast.InternalGains(:,i) = b_loads.InternalGains;
        b_forecast.NonHVACelectric(:,i) = b_loads.Equipment + b_loads.InteriorLighting + b_loads.ExteriorLighting + b_loads.OtherLoads;
        if isfield(b_loads,'DCloads')
            b_forecast.DCloads = b_loads.DCloads;
        end
    end
end

%%Need warm-up period if not currently running the model
for i = 1:1:n_b
    if  abs(round(864000*(buildings(i).Timestamp+options.Resolution/24))/864000 - date(1))>1e-5
        if length(date) == 1
            wu_date = linspace(date(1) + options.Resolution/24,date(1)+1,24)';
        else
            wu_date = linspace(date(1),date(1)+1 - options.Resolution/24,24)';
        end
        if isempty(hist_prof)
            f_names = fieldnames(prev_data.Weather);
            for f = 1:1:length(f_names)
                wu_weather.(f_names{f}) = interp1(prev_data.Timestamp,prev_data.Weather.(f_names{f}),wu_date);
            end
        else
            wu_weather = weather_forecast(prev_data,hist_prof,wu_date);
        end
        buildings(i) = building_warmup(buildings(i),wu_weather,options.Resolution/24,wu_date,6);
    end
end
b_forecast.E0 = zeros(n_s,n_b);
b_forecast.C0 = zeros(n_s,n_b);
b_forecast.H0 = zeros(n_s,n_b);
b_forecast.Tmin = zeros(n_s,n_b);
b_forecast.Tmax = zeros(n_s,n_b);
b_forecast.Tset_H = zeros(n_s,n_b);
b_forecast.Tset_C = zeros(n_s,n_b);
b_forecast.Fan_Power = zeros(n_s,n_b);
b_forecast.AirFlow = zeros(n_s,n_b);
b_forecast.Tzone = zeros(n_s,n_b);
b_forecast.Twall = zeros(n_s,n_b);
b_forecast.Damper = zeros(n_s,n_b);
for i = 1:1:n_b
    [cooling, heating, fan_power,zone,wall,damper] = building_profile(buildings(i),date,b_forecast.InternalGains(:,i),b_forecast.ExternalGains(:,i),weather.Tdb,weather.RH,buildings(i).Tzone,buildings(i).Twall);
    b_forecast.Tmin(:,i) = load_sched(buildings(i),date,'TsetH') - 0.5*buildings(i).VariableStruct.Comfort;
    b_forecast.Tmax(:,i) = load_sched(buildings(i),date,'TsetC') + 0.5*buildings(i).VariableStruct.Comfort;
    [temperature_to_heat,temperature_to_cool] = equilib_temperature(buildings(i),date,weather,zone(2:end),wall(2:end),cooling,heating,b_forecast.InternalGains(:,i),b_forecast.ExternalGains(:,i),b_forecast.Tmin(:,i),b_forecast.Tmax(:,i));
    b_forecast.E0(:,i) = b_forecast.NonHVACelectric(:,i) + fan_power(:,i);
    b_forecast.C0(:,i) = cooling(:,i) ;
    b_forecast.H0(:,i) = heating(:,i);
    if ~buildings(i).QPform.Heating
        %Can upgrade to include non-linear mapping of heating to electricity for built-in HVAC equipment
        b_forecast.E0(:,i) = b_forecast.E0(:,i) + (b_forecast.H0(:,i)- buildings(i).QPform.UA*(zone(2:end)-temperature_to_heat))*buildings(i).QPform.H2E;
    end
    if ~buildings(i).QPform.Cooling
        %Can upgrade to include non-linear mapping of cooling to electricity for built-in HVAC equipment
        b_forecast.E0(:,i) = b_forecast.E0(:,i) + (b_forecast.C0(:,i)- buildings(i).QPform.UA*(temperature_to_cool-zone(2:end)))*buildings(i).QPform.C2E;
    end
    b_forecast.Tset_H(:,i) = temperature_to_heat;
    b_forecast.Tset_C(:,i) = temperature_to_cool;
    b_forecast.Fan_Power(:,i) = fan_power;
    b_forecast.AirFlow(:,i) = fan_power/buildings(i).VariableStruct.FanPower;
    b_forecast.Tzone(:,i) = zone(2:end);
    b_forecast.Twall(:,i) = wall(2:end);
    b_forecast.Damper(:,i) = damper;
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
dry_bulb_kelvin = weather.Tdb+273.15; %Tdb (Kelvin)
sat_press = exp((-5.8002206e3)./dry_bulb_kelvin + 1.3914993 - 4.8640239e-2*dry_bulb_kelvin + 4.1764768e-5*dry_bulb_kelvin.^2 - 1.4452093e-8*dry_bulb_kelvin.^3 + 6.5459673*log(dry_bulb_kelvin))/1000; %saturated water vapor pressure ASHRAE 2013 fundamentals eq. 6 in kPa valid for 0 to 200C
water_press = weather.RH/100.*sat_press; % kPa
amb_spec_heat = 1.006 + 1.86*(.621945*(water_press./(pressure-water_press))); % kJ/kg*K
min_flow = build.VariableStruct.Volume*rho_Air*build.VariableStruct.AirChangePerHr*(1/3600); %kg/s of flow
for t = 1:1:n_s
    temperature_range = linspace(min_temp(t)-2*build.VariableStruct.Comfort,max_temp(t)+2*build.VariableStruct.Comfort,n)';
    net_gain = zone_thermal_load(build,wall(t),weather.Tdb(t),temperature_range,internal_gains(t),external_gains(t),dt(t));
    need_heat = (min(temperature_range,build.VariableStruct.ColdAirSet) - weather.Tdb(t))*amb_spec_heat(t)*min_flow*dt(t) - net_gain;
    if all(need_heat>0)
        temperature_to_heat(t) = temperature_range(1);
    elseif all(need_heat<0)
        temperature_to_heat(t) = temperature_range(end);
    else
        temperature_to_heat(t) = interp1(need_heat,temperature_range,0);
    end
    temperature_to_heat = max(temperature_to_heat,zone(t)-heating(t)*build.QPform.UA);
    need_cool = (weather.Tdb(t) - temperature_range)*amb_spec_heat(t)*min_flow*dt(t) + net_gain;
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