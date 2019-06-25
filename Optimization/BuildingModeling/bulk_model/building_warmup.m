function buildings = building_warmup(buildings,options,date,prev_data,hist_prof)
%%Need warm-up period if not currently running the model
n_b = length(buildings);
wu_date = linspace(date(1) + options.Resolution/24,date(1)+1,24)';
if isempty(hist_prof)
    f_names = fieldnames(prev_data.Weather);
    for f = 1:1:length(f_names)
        wu_weather.(f_names{f}) = interp1(prev_data.Timestamp,prev_data.Weather.(f_names{f}),wu_date);
    end
else
    wu_weather = weather_forecast(prev_data,hist_prof,wu_date);
end
for i = 1:1:n_b
    if  ~isfield(buildings(i),'Tzone') || isempty(buildings(i).Tzone) || abs(round(864000*(buildings(i).Timestamp+options.Resolution/24))/864000 - date(1))>1e-5
        zone = 20;
        wall = 20;
        if length(date) == 1
            wuDate = linspace(date(1) + options.Resolution/24,date(1)+1,24)';
        else
            wuDate = linspace(date(1),date(1)+1 - options.Resolution/24,24)';
        end
        wu_sg = solar_gain(buildings(i),wuDate,buildings(i).QPform.Location,wu_weather);
        wu_loads = building_loads(buildings(i),wuDate,wu_sg);
        wu_external_gains = wu_sg.Walls + wu_sg.Roof;
        for d = 1:1:6
            [~,~,~,zone,wall,~] = building_profile(buildings(i),wuDate,wu_loads.InternalGains,wu_external_gains,wu_weather.DrybulbC,wu_weather.RHum,zone(end,1),wall(end,1));
        end
        if zone(end)<10 || zone(end)>40
            disp('BuildingWarmUp Error')
        end
        buildings(i).Tzone = zone(end);
        buildings(i).Twall = wall(end);
        buildings(i).Timestamp = date(1);
    end
end
end%Ends function building_warmup