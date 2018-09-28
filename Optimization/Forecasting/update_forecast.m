function [forecast,gen,buildings] = update_forecast(gen,buildings,fluid_loop,subnet,options,date,test_data,dispatch)
freq = max(1,round(options.Horizon/24)); %period of repetition (1 = 1 day) This is how far back the forecasting methods are able to see. It is irrelevant if the forecast is perfect
res = options.Resolution/24;
n_o = round(freq/res)+1;
prev_data = get_data(test_data,linspace((date(1) - 2*res - freq),date(1)-2*res,n_o)',[],[]);

if isfield(test_data,'HistProf')
    hist_prof = test_data.HistProf;
else
    hist_prof = [];
end
if length(date) == 1
    forecast = get_data(test_data,date(1),[],[]);
    if ~isempty(buildings)
        buildings = building_warmup(buildings,options,date,prev_data,hist_prof);
    end
else
    switch options.forecast
        case 'ARMA'
            forecast = arma(date,prev_data);
            forecast.Weather = weather_forecast(prev_data,hist_prof,date);
        case 'ARIMA'
            forecast = arima_wsu(date,prev_data,options);
            forecast.Weather = weather_forecast(prev_data,hist_prof,date);
        case 'NeuralNet'
            %
        case 'Surface'
            weather = weather_forecast(prev_data,hist_prof,date);
            forecast = surface_forecast(prev_data,hist_prof.Demand,date,weather.DrybulbC,[]);
            forecast.Weather = weather;
        case 'Perfect'
            forecast = get_data(test_data,date,[],[]);
        case 'Building'
            forecast.Timestamp = date;
            forecast.Weather = weather_forecast(prev_data,hist_prof,date);
    end
end
if strcmp(options.method,'Dispatch') && length(date)>1
    %make 1st hour forecast "perfect"
    next_data = get_data(test_data,date(1),[],[]);
    f = fieldnames(next_data);
    for j = 1:1:length(f)
        if isstruct(next_data.(f{j}))
            S = fieldnames(next_data.(f{j}));
            for i = 1:1:length(S)
                forecast.(f{j}).(S{i})(1,:) = next_data.(f{j}).(S{i})(1,:);
            end
        else
            forecast.(f{j})(1,:) = next_data.(f{j})(1,:);
        end
    end
end
if ~isempty(buildings) 
    [buildings,forecast] = forecast_building(forecast,date,buildings,options);
end
if isfield(forecast,'Weather') && isfield(forecast.Weather,'DNIWm2')
    forecast.Renewable = renewable_output(gen,subnet,date,forecast.Weather.DNIWm2);
end
if isfield(subnet,'Hydro')
    now_data = get_data(test_data,date(1)-res,[],[]);
    if ~isempty(dispatch) && isfield(dispatch,'LineFlows')
        xi = max(1,nnz((dispatch.Timestamp-1e-6)<prev_data.Timestamp(1) & dispatch.Timestamp>0));
        xf = nnz((dispatch.Timestamp-1e-6)<=now_data.Timestamp & dispatch.Timestamp>0);
        yi = 1 + nnz(prev_data.Timestamp<dispatch.Timestamp(1));
        prev_data.Hydro.OutFlow = zeros(size(prev_data.Hydro.InFlow));
        now_data.Hydro.OutFlow = zeros(size(now_data.Hydro.InFlow));
        for n = 1:1:length(subnet.Hydro.lineNumber)
            prev_data.Hydro.OutFlow(yi:end,n) = interp1(dispatch.Timestamp(xi:xf),dispatch.LineFlows(xi:xf,subnet.Hydro.lineNumber(n)),prev_data.Timestamp(yi:end));
            now_data.Hydro.OutFlow(n) = dispatch.LineFlows(xf,subnet.Hydro.lineNumber(n));
        end
        forecast.Hydro.InFlow = forecast_hydro(date,forecast.Hydro.SourceSink,subnet,prev_data,now_data);
    else
        %%estimate history for InFlow
    end
end

if options.SpinReserve
    if isfield(forecast,'Demand')
        forecast.SRtarget = forecast.SRtarget + options.SpinReservePerc/100*sum(forecast.Demand.E,2);
    end
end
end%Ends function update_forecast

function in_flow = forecast_hydro(date,source_sink,subnet,prev_data,now_data)
% forecast the source/sink + previously released upsteam water flow at each 
% time step: 1 to nS. When t<T it needs the scheduled outflow of the
% upstream dam to be on the constant side of the optimization.
% does not give a forecast at t=0;
n_s = length(date);
t_last_disp = now_data.Timestamp;
in_flow = source_sink;
time = [prev_data.Timestamp;now_data.Timestamp];
outflow = [prev_data.Hydro.OutFlow;now_data.Hydro.OutFlow];
for n = 1:1:length(subnet.Hydro.nodes) 
    K = subnet.Hydro.up_river{n};
    for j = 1:1:length(K)
        T_seg = subnet.Hydro.lineTime(K(j));
        for t = 1:1:n_s
            if (date(t) - t_last_disp)<= T_seg/24 % All of the flow reaching this reservior at this time was previously scheduled
                in_flow(t,n) = in_flow(t,n) + interp1(time,outflow(:,K(j)),date(t)-T_seg/24); %outflow from previous dispatches because the time preceeds t_last_disp
            elseif t>1 && (date(t-1) - t_last_disp)< T_seg/24 % A portion of the flow reaching this reservior at this time was previously scheduled
                frac = (date(t)- t_last_disp - T_seg/24)/(date(t) - date(t-1));
                if (date(t)-T_seg/24)>time(end)
                    in_flow(t,n) = in_flow(t,n) + (1-frac)*outflow(end,K(j));
                else
                    in_flow(t,n) = in_flow(t,n) + (1-frac)*interp1(time,outflow(:,K(j)),date(t)-T_seg/24); %outflow from previous dispatches because the time preceeds t_last_disp
                end
            end
            if isnan(in_flow(t,n))
                disp('Potential rounding error in update_forecast')
            end
        end
    end
end
end %Ends function forecast_hydro