function forecast = surface_forecast(prev,hist_prof,date,temperature_forecast,holidays)
forecast.Timestamp = date;
d_sim = prev.Timestamp(end);
days = ceil(date(end)-floor(d_sim)); %this says that you may need to use surface fits for multiple days.
s12hour = nnz(date<(d_sim+.5)); %steps in next 12 hours
a = datevec(d_sim);
date_days = floor(d_sim)-1:ceil(d_sim+days);
mon_fri = zeros(length(date_days),1);
for i = 1:1:length(date_days)
    if weekday(date_days(i))<=6 && weekday(date_days(i))>=2 && (isempty(holidays) || ~any(holidays==date_days(i)))
        mon_fri(i) = 1;
    else
        mon_fri(i) = 0;
    end
end
month_name = {'Jan' 'Feb' 'Mar' 'Apr' 'May' 'Jun' 'Jul' 'Aug' 'Sep' 'Oct' 'Nov' 'Dec' 'Jan'};
hour = (date-floor(d_sim))*24;
f = fieldnames(prev.Demand);
for s = 1:1:length(f) %repeat for electric, cooling, heating, and steam as necessary
    [~,n_d] = size(prev.Demand.(f{s}));
    forecast.Demand.(f{s}) = zeros(length(date),n_d);
    list = fieldnames(hist_prof.(f{s})(1));
    for k = 1:1:n_d
        month = a(2);
        hist_fit_dem = nan(length(date),1);
        for j = 0:1:ceil(days)
            if a(3)+ (j-1) > (datenum(a(1), a(2)+1, 1)-datenum(a(1), a(2), 1))
                month = a(2)+1;
            end
            %load suface
            if length(list)>1
                if nnz(strcmp(list,strcat(month_name(month),'WeekEnd')))>0 %historical profile is split to weekday/weekend
                    if mon_fri(j+1) == 0
                        surface = hist_prof.(f{s})(k).(strcat(month_name{month},'WeekEnd'));
                    else
                        surface = hist_prof.(f{s})(k).(strcat(month_name{month},'WeekDay'));
                    end
                elseif nnz(strcmp(list,char(month_name(month))))>0 %historical profile is broken by month, but not weekend/weekday
                    surface = hist_prof.(f{s})(k).(month_name{month});
                else 
                    surface = hist_prof.(f{s})(k).(list{1});
                end
            else
                surface = hist_prof.(f{s})(k).(list{1});% only 1 historical surface fit
            end
            if j ==0 %surface for yesterday
                yest_pred = surface(mod((prev.Timestamp-floor(prev.Timestamp(1)))*24,24),prev.Weather.Tdb);
                error_pred = min(1,max(-1,(prev.Demand.(f{s})(:,k)- yest_pred)./yest_pred));
                bias_perc = mean(error_pred)*100; %average percent error between hisorical expectations for yesterday, and actual yesterday
                hourly_error = prev.Demand.(f{s})(:,k) - yest_pred*(1+0.7*bias_perc/100);%hourly variation with equal mean (historical and last 24 hours)
            else %surface for current day or subsequent days
                index = (hour<=(24*j)) & (hour>=(24*(j-1))); %index of Time corresponding to this day
                hist_fit_dem(index) = surface(mod(hour(index),24),temperature_forecast(index));
            end
        end
        forecast.Demand.(f{s})(1:s12hour,k) = max(0,hist_fit_dem(1:s12hour)*(1+0.7*bias_perc/100) + linspace(1,0,s12hour)'.*interp1(prev.Timestamp,hourly_error,date(1:s12hour)-1));
        forecast.Demand.(f{s})(s12hour+1:end,k) = max(0,hist_fit_dem(s12hour+1:end)*(1+0.7*bias_perc/100));
        forecast.Demand.(f{s})(abs(forecast.Demand.(f{s}))<1e-2) = 0; 
    end
end
end%ends function surface_forecasts