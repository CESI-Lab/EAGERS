function solar = load_solar(timestamp)
dir = strrep(which('load_solar.m'),fullfile('Optimization','SetupFunctions','load_solar.m'),'');
cd(fullfile(dir,'Data','Solar'))
[fn,pn,~] = uigetfile('*.mat','Load Solar Irradiation Data File');
load(fullfile(pn,fn));
solar = zeros(length(timestamp),1);
if timestamp(1)>=weather.Timestamp(1) && timestamp(end)<=weather.Timestamp(end)
    solar = interp1(weather.Timestamp,weather.irradDireNorm,timestamp); 
else %need to move solar data timestamp to line up with timestamp
    D1 = datevec(timestamp(1));
    D2 = datevec(weather.Timestamp(1));
    if datenum([D1(1),D2(2),D2(3),D2(4),D2(5),D2(6)])>timestamp(1)
        D2(1) = D2(1)-1;
    end
    weather.Timestamp = weather.Timestamp + datenum([D1(1),1,1]) - datenum([D2(1),1,1]);%days between start of years for testdata and solar data respectively
    y = 1;
    Xi = 1;
    Xf = nnz(timestamp<=weather.Timestamp(end));
    while Xi<=length(timestamp)
        solar(Xi:Xf,1) = interp1(weather.Timestamp,weather.irradDireNorm,timestamp(Xi:Xf)); 
        weather.Timestamp = weather.Timestamp + datenum([D1(1)+y,1,1]) - datenum([D1(1)+y-1,1,1]);%shift weather data 1 year
        Xi = Xf+1;
        Xf = nnz(timestamp<=weather.Timestamp(end));
        y = y+1;
    end
end
cd(dir)

end%ends function load_solar