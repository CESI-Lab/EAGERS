%click on any weather file in the weather folder
global  TestData
options.tspacing='constant';
options.Horizon = 24;
options.Resolution = 1; 
nS = options.Horizon/options.Resolution;
Time = linspace(options.Resolution,options.Horizon,nS)';
Date = linspace(datenum([2017,1,1]),datenum([2018,1,1]),8760/options.Resolution+1)';

TestData.Timestamp = Date;
TestData.Weather = interpolate_weather(weather,Date);
S = fieldnames(TestData.Weather);
for j = 1:1:length(S)
    if isnumeric(TestData.Weather.(S{j}))
        HistProf.(S{j}) = typical_day([],Date,TestData.Weather.(S{j}));
        Forecast.(S{j}) = zeros(length(Date),round(24/options.Resolution));
    end
end
for t = 1:1:length(Date)
    forecasttime = Date(t) + Time;
    W = weather_forecast(TestData,HistProf,forecasttime,options.Resolution);
    for j = 1:1:length(S)
        if isnumeric(TestData.Weather.(S{j}))
            Forecast.(S{j})(t,:) = W.(S{j})';
        end
    end
end