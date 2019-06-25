function hist_prof = calculate_fits(data,timestamp,weather)
f = fieldnames(data);
for s = 1:1:length(f)
    dem = f{s};
    h=waitbar(0,'Recalculating surface fit');

    steps = round(1/(timestamp(2)-timestamp(1)));%points in 1 day of data
    resolution =24/steps;
    [~,m,~,hour,minutes,~] = datevec(timestamp(2:end-1));%avoid checking month 1st and last point if starting or ending at 00
    m = [m(1); m; m(end)];%make m the same length as data
    hour = [max(0,floor(hour(1)-resolution)); hour; min(24,hour(end)+floor((minutes(end)+60*resolution)/60))];%make hour the same length as data
    minutes = [max(0,minutes(1)-60*resolution); minutes; mod(minutes(end)+60*resolution,60)];%make the same length as data
    months = sort(unique(m));

    monthNames = cellstr({'Jan';'Feb';'Mar';'Apr';'May';'Jun';'Jul';'Aug';'Sep';'Oct';'Nov';'Dec';});
    bus_day = true(length(timestamp),1);
    days_after1_1_17 = floor(timestamp) - datenum([2017,1,1]);
    wd = 1 + mod(days_after1_1_17,7); %day of the week, sunday is 1, saturday is 7
    bus_day(wd==1) = false;
    bus_day(wd==7) = false;

    if length(months)<=4 %make 1 surface fit
        for k = 1:1:length(data.(dem)(1,:))
            z = data.(dem)(:,k);
            y = weather.Tdb;
            x = hour+minutes/60;
            if steps>100 % too many data points for a surface fit, average each hour
                r = ceil(steps/24);
                n2 = floor(length(x)/r);
                x2 = zeros(n2,1);
                y2 = zeros(n2,1);
                z2 = zeros(n2,1);
                for i=1:1:n2%turn while loop into for loop to make loop faster
                    x2(i)=sum(x((i-1)*r+1:i*r))/r;
                    y2(i)=sum(y((i-1)*r+1:i*r))/r;
                    z2(i)=sum(z((i-1)*r+1:i*r))/r;
                end
                x = x2;
                y = y2;
                z = z2;
            end
            hist_prof.(dem)(k).Annual = fit([x y],z,'lowess');
        end
    else %split surface fit by month and weekday/weekend
        y = weather.Tdb;
        x = hour+minutes/60;
        for k = 1:1:length(data.(dem)(1,:))
            z = data.(dem)(:,k);
            for i = 1:1:length(months)
                waitbar((i-1)/length(months),h,strcat('Calculating surface fit for ',dem,' for the month of ',monthNames{months(i)}));
                if steps>150 % too many data points for a surface fit, average each 10 min
                    r = ceil(steps/(24*6));
                    n2 = floor(length(x)/r);
                    x2=zeros(n2,1);
                    y2=zeros(n2,1);
                    z2=zeros(n2,1);
                    for j=1:1:n2%turn the while loop into a for loop to make loop faster
                        x2(j)=sum(x((j-1)*r+1:j*r))/r;
                        y2(j)=sum(y((j-1)*r+1:j*r))/r;
                        z2(j)=sum(z((j-1)*r+1:j*r))/r;
                    end
                    x = x2;
                    y = y2;
                    z = z2;
                end
                % business days
                h1 = (m==months(i)) & bus_day; %points in the correct month
                z1 = z(h1);%remove the days that are non business days
                y1 = y(h1);
                x1 = x(h1);%cannot use nonzeros function because there are some elements in X1 that are 0
                % weekends/holidays
                h2 = (m==months(i)) & ~bus_day;
                z2 = z(h2);
                y2 = y(h2);
                x2 = x(h2);
                hist_prof.(dem)(k).(strcat(monthNames{months(i)},'WeekDay')) = fit([x1, y1],z1,'lowess');%'poly23');%'loess');
                hist_prof.(dem)(k).(strcat(monthNames{months(i)},'WeekEnd')) = fit([x2, y2],z2,'lowess');%'poly23');%'loess');
            end
        end
    end
    close(h)
end
end%ends function calculate_fits