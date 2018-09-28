function test_weather = interpolate_weather(weather,date)
test_weather.filename = weather.filename;
if ~isfield(weather,'Timestamp')%Assume it is a typical meteorological year
    weather.Timestamp = linspace(datenum([2017,1,1,8760/length(weather.DrybulbC),0,0]),datenum([2018,1,1]),length(weather.DrybulbC))';
end
n_s = length(date);
s = fieldnames(weather);
s = s(~strcmp('Timestamp',s));
for j = 1:1:length(s)
    if isnumeric(weather.(s{j}))
        test_weather.(s{j}) = nan(n_s,1);
    end
end
if date(1)>=weather.Timestamp(1) && date(end)<=weather.Timestamp(end)
    for j = 1:1:length(s)
        if isnumeric(weather.(s{j}))
            test_weather.(s{j}) = interp1(weather.Timestamp,weather.(s{j}),date); 
        end
    end
else %need to move solar data timestamp to line up with TestData timestamp
    d1 = datevec(date(1));
    d2 = datevec(weather.Timestamp(1));
    
    weather.Timestamp = weather.Timestamp + datenum([d1(1),1,1]) - datenum([d2(1),1,1]);%days between start of years for testdata and solar data respectively
    if datenum([d1(1),d2(2),d2(3),d2(4),d2(5),d2(6)])>date(1)
        weather.Timestamp = weather.Timestamp - (datenum([d1(1),1,1]) - datenum([d1(1)-1,1,1]));
    end
    y = 1;
    while y == 1 || weather.Timestamp(1)<=date(end)
        index = date<=weather.Timestamp(end) & date>=weather.Timestamp(1);
        for j = 1:1:length(s)
            if isnumeric(weather.(s{j}))
                test_weather.(s{j})(index,1) = interp1(weather.Timestamp,weather.(s{j}),date(index)); 
            end
        end
        d3 = datevec(weather.Timestamp(1));
        weather.Timestamp = weather.Timestamp + datenum([d3(1)+1,1,1]) - datenum([d3(1),1,1]);%shift weather data 1 year
        y = y+1;
    end
    for j = 1:1:length(s)
        if isnumeric(weather.(s{j}))
            index = isnan(test_weather.(s{j}));
            test_weather.(s{j})(index,1) = weather.(s{j})(1);
        end
    end
end