function forecast = weather_forecast(prev,hist_prof,date)
forecast = [];
if isfield(prev,'Weather')
    A = datevec(prev.Timestamp(end));
    hour = (date-floor(prev.Timestamp(end)))*24;
    F = fieldnames(prev.Weather);
    for j = 1:1:length(F) 
        if isnumeric(prev.Weather.(F{j}))
            if length(date) == 1
                forecast.(F{j}) = prev.Weather.(F{j})(1);
            else
                YestFit = interp1(linspace(0,24,length(prev.Weather.(F{j})))',prev.Weather.(F{j}),mod((date-prev.Timestamp(end))*24,24));
                HistFit = interp1(0:24,[hist_prof.(F{j})(A(2),end),hist_prof.(F{j})(A(2),:)],mod(hour,24));
                W = interp1([prev.Timestamp(end),prev.Timestamp(end)+3,prev.Timestamp(end)+100],[0.9,0,0],date);%weight between yesterday and historical average
                W(isnan(W)) = 0;
                forecast.(F{j}) =  (W.*YestFit + (1-W).*HistFit);% Balanced between yesterdays T and historical T (includes forecast for current time)
            end
        end
    end
end
end %ends weather_forecast