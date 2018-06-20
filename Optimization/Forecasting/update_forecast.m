function [forecast,gen,buildings] = update_forecast(gen,buildings,cool_tower,subnet,options,date,hist_prof,prev_data,now_data,future_data)
%Date is the date number, Time is a vector of times (in hours)
if length(date) == 1
    forecast = now_data;
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
            Weather = weather_forecast(prev_data,hist_prof,date);
            forecast = surface_forecast(prev_data,hist_prof.Demand,date,Weather.Tdb,[]);
            forecast.Weather = Weather;
        case 'Perfect'
            forecast = future_data;
        case 'Building'
            forecast.Timestamp = date;
            forecast.Weather = weather_forecast(prev_data,hist_prof,date);
    end
end
if strcmp(options.method,'Dispatch') && length(date)>1
    %make 1st hour forecast "perfect"
    f = fieldnames(future_data);
    for j = 1:1:length(f)
        if isstruct(future_data.(f{j}))
            S = fieldnames(future_data.(f{j}));
            for i = 1:1:length(S)
                forecast.(f{j}).(S{i})(1,:) = future_data.(f{j}).(S{i})(1,:);
            end
        else
            forecast.(f{j})(1,:) = future_data.(f{j})(1,:);
        end
    end
end
if ~isempty(buildings) && ~strcmp(options.solver,'NREL')
    [buildings,forecast.Building] = forecast_building(forecast.Weather,date,buildings,prev_data,hist_prof,options);
end
if isfield(forecast,'Weather') && isfield(forecast.Weather,'irradDireNorm')
    forecast.Renewable = renewable_output(gen,subnet,date,forecast.Weather.irradDireNorm);
end
if isfield(subnet,'Hydro')
    forecast.Hydro.InFlow = forecast_hydro(date,forecast.Hydro.SourceSink,subnet,prev_data,now_data);
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