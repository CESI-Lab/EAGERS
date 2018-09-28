function profile = load_sched(schedules,holidays,date,param,dd)
sd = floor(date(1)); % start date
n_s = length(date); % number timesteps
daylight_savings = find_daylight_savings(date(1));
holiday = find_holidays(date,holidays);
hour = round(1000*24*(date-floor(date)))/1000;
hour(date==floor(date)) = 24; % change hour 0 to hour 24 for interpolating.
daysAfter1_1_17 = floor(date(1)) - datenum([2017,1,1]);
if isempty(param)
    param = fieldnames(schedules);
    single_sched = false;
else
    single_sched = true;
end
for i = 1:1:length(param)
    sched = schedules.(param{i});
    profile.(param{i}) = zeros(n_s,1);
    if isfield(sched,'Timestamp')
        profile.(param{i}) = interp_sched(sched,date);%user defined schedules defined at each time step
    else
        wd = 1 + mod(daysAfter1_1_17,7); % day of the week, sunday is 1, saturday is 7
        p2 = 0;
        ds = 0;
        while p2<n_s
            p = p2+1;
            d = floor(date(p)+ds/24-date(1))+1;
            if (d>=daylight_savings.start && d<=daylight_savings.end) || (daylight_savings.end<daylight_savings.start && (d<=daylight_savings.end || d>daylight_savings.start))
                ds = 1;
            else
                ds = 0;
            end
            p2 = p;
            while p2<n_s && date(p2+1)<=sd+d-ds/24
                p2 = p2+1;
            end
            if ~isempty(dd) && strcmp(dd,'summer') && isfield(sched,'SummerDesignDay')
                day = 'SummerDesignDay';
            elseif ~isempty(dd) && strcmp(dd,'winter') && isfield(sched,'WinterDesignDay')
                day = 'WinterDesignDay';
            elseif isempty(dd) && holiday(d) && isfield(sched,'Holidays')
                day = 'Holidays';
            elseif isempty(dd) && holiday(d) && isfield(sched,'AllOtherDays')
                day = 'AllOtherDays';
            elseif isfield(sched,'AllDays')
                day = 'AllDays';
            elseif isfield(sched,'ALLDays')
                day = 'ALLDays';
            else
                if wd == 1 && isfield(sched,'Sunday')% Sunday
                    day = 'Sunday';
                elseif wd == 7 && isfield(sched,'Saturday')% Saturday
                    day = 'Saturday';
                    wd = wd-7;
                elseif wd>1 && wd < 7 
                    if isfield(sched,'Weekday')% Weekday
                        day = 'Weekday';
                    elseif isfield(sched,'Weekdays')% Weekdays
                        day = 'Weekdays';
                    else
                        day = 'AllOtherDays';
                    end
                else
                    day = 'AllOtherDays';
                end
            end
            if length(sched.seasons)>1
                %correct for season
                d_now = datevec(date(p));
                h_of_y = 24*(date(p) - datenum([d_now(1),1,1])); % hour since start of year
                s = nnz(sched.seasons<=h_of_y)+1;% get season
            else
                s=1;
            end
            interp_hour = hour(p:p2) + ds;
            interp_hour(interp_hour>24) = interp_hour(interp_hour>24)-24;
            if sched.interpolate(s) && sched.ramp<1e-3 %EnergyPlus type interpolation
                profile.(param{i})(p:p2,1) = interp1(sched.(day){s}(:,1),sched.(day){s}(:,2),interp_hour);
            else
                sched_2_interp = convert_sched(sched.(day){s},sched.ramp);
                profile.(param{i})(p:p2,1) = interp1(sched_2_interp(:,1),sched_2_interp(:,2),interp_hour);
            end
            wd = wd+1;
        end
    end
end
if single_sched
    profile = profile.(param);
end
end%Ends function load_sched

function profile = interp_sched(sched,date)
%Interpolates two fields Data.Timestamp and Data.Load
%If the date is in the timestamp range, it uses the value
%If not it tries to use the same day of the 1st year in the timestamp
%If that day doesn't exist for a different year, it maches the day of the month,

days = max(1,ceil(date(end) - floor(date(1))));
n_s = length(date); % number timesteps
profile = zeros(n_s,1);

dt = date(2)-date(1);
[newDay,D_of_Y, D_of_M] = day_of_year(date);
r = (sched.Timestamp(2)-sched.Timestamp(1))/dt;% # of points that must be created between datum 
for d = 1:1:days
    p = newDay(d);
    if d<days
        p2 = newDay(d+1)-1;
    else 
        p2 = length(date);
    end
    if date(p)>=floor(sched.Timestamp(1)) && (date(p2)-dt)<floor(sched.Timestamp(end))
        %Date is within the data period, use directly
        if sched.Timestamp(1)>(date(p)-dt)
            Xi = 1;
        else
            Xi = nnz(sched.Timestamp<=date(p)-dt);
        end
        Xf = nnz(sched.Timestamp<=(date(p)+1));
    else%find matching day in different year
        [data_newDay,data_D_of_Y,data_D_of_M] = day_of_year(sched.Timestamp);
        data_days = length(data_newDay);
        if any(data_D_of_Y == D_of_Y(d))
            Index = nonzeros((1:data_days)'.*(data_D_of_Y == D_of_Y(d)));
            Xi = data_newDay(Index(1));
            if Index(1) == data_days
                Xf = length(sched.Timestamp);
            else
                Xf = data_newDay(Index(1)+1)-1;
            end
        elseif any(data_D_of_M == D_of_M(d))
            Index = nonzeros((1:data_days)'.*(data_D_of_M == D_of_M(d)));
            Xi = data_newDay(Index(1));
            if Index(1) == data_days
                Xf = length(sched.Timestamp);
            else
                Xf = data_newDay(Index(1)+1)-1;
            end
        else
            rn = ceil(rand(1)*data_days);
            Xi = data_newDay(rn);
            if rn == data_days
                Xf = length(sched.Timestamp);
            else
                Xf = data_newDay(rn+1)-1;
            end
        end
    end
    if abs(r-1)<1e-8
        profile(p:p2,1) = sched.Load(Xi:Xf);
    elseif r<1 %extra datum, average points in between
        nS2 = length(sched.Timestamp);
        x1 = Xi;
        for k = p:1:p2
            if k == p2
                x2 = Xf;
            else
                x2 = x1;
                while x2<nS2 && (sched.Timestamp(x2+1)-floor(sched.Timestamp(x2+1)))<(date(k)-floor(date(k)))
                    x2 = x2+1;
                end
            end
            profile(k,1) = mean(sched.Load(x1:x2));
            x1 = x2+1;
        end
    elseif r>1 %interpolate between timesteps
        t = sched.Timestamp(Xi:Xf) - floor(sched.Timestamp(Xi:Xf));
        t2 = date(p:p2)-floor(date(p));
        profile(p:p2,1) = interp1(t,sched.Load(Xi:Xf),t2);
    end
end
end%Ends function interp_sched

function [newDay,D_of_Y, D_of_M] = day_of_year(Date)
D = datevec(Date(1));
sd = floor(Date(1)); % start date
days = max(1,ceil(Date(end) - datenum([D(1),D(2),D(3)])));
newDay = zeros(days,1);
D_of_Y = zeros(days,1);
D_of_M = zeros(days,1);
p = 1;
nS = length(Date);
for d = 1:1:days
    p2 = p;
    while p2<nS && Date(p2+1)<=sd+d
        p2 = p2+1;
    end
    newDay(d) = p;
    D = datevec(Date(p));
    D_of_Y(d) = floor(Date(p) - datenum([D(1),1,1]))+1;
    D_of_M(d) = floor(Date(p) - datenum([D(1),D(2),1]))+1;
    p = p2+1;
end
end%Ends function day_of_year

function daylight_savings = find_daylight_savings(date)
[Y,~,~] = datevec(floor(date));
yd = datenum([Y,1,1]) + (0:365)';
D = datevec(yd);
wd = 1 + mod(yd - datenum([2017,1,1]),7); % day of the week, sunday is 1, saturday is 7
march_sunday = nonzeros((1:366)'.*(D(:,2)==3 & wd == 1));
daylight_savings.start = march_sunday(2) - (floor(date) - datenum([Y,1,1]));
november_sunday = nonzeros((1:366)'.*(D(:,2)==11 & wd == 1));
daylight_savings.end = november_sunday(1) - (floor(date) - datenum([Y,1,1]));
end%Ends function find_daylight_savings

function holiday = find_holidays(date,holidays)
days = ceil(date(end) - floor(date(1)))+1;
date_days = (floor(date(1)):floor(date(1))+days-1)';
D = datevec(date_days);
jan1 = datenum([D(1,1),1,1]);
d_o_y = mod(date_days - jan1+1,365);
daysAfter1_1_17 = (date_days - datenum([2017,1,1]));
wd = 1 + mod(daysAfter1_1_17,7); % day of the week, sunday is 1, saturday is 7
holiday = false(days,1);
month_names = {'January';'February';'March';'April';'May';'June';'July';'August';'September';'October';'November';'December';};
month_days = [31;28;31;30;31;30;31;31;30;31;30;31];
weekday_names = {'Sunday';'Monday';'Tuesday';'Wednesday';'Thursday';'Friday';'Saturday';};
for i = 1:1:length(holidays.name)
    spc = strfind(holidays.start{i},' ');
    seg1 = holidays.start{i}(1:spc(1)-1);
    month = nonzeros((1:12)'.*strcmp(seg1,month_names));
    if ~isempty(month)
        day = str2double(holidays.start{i}(spc(1)+1:end));
        h_d = sum(month_days(1:month-1)) + day;
        holiday(d_o_y==h_d) = true;
    else
        d_o_w = nonzeros((1:7)'.*strcmp(holidays.start{i}(spc(1)+1:spc(2)-1),weekday_names));
        seg4 = holidays.start{i}(spc(3)+1:end);
        month = nonzeros((1:12)'.*strcmp(seg4,month_names));
        switch holidays.start{i}(1:spc(1)-1)
            case 'Last'
                pos_days = nonzeros((1:days)'.*(D(:,2) == month & wd== d_o_w & D(:,3)>21));
                if ~isempty(pos_days)
                    holiday(pos_days(end)) = true;
                end
            case '1st'
                pos_days = nonzeros((1:days)'.*(D(:,2) == month & wd== d_o_w & D(:,3)<=7));
                holiday(pos_days) = true;
            case '2nd'
                pos_days = nonzeros((1:days)'.*(D(:,2) == month & wd== d_o_w & D(:,3)>7 & D(:,3)<=14));
                holiday(pos_days) = true;
            case '3rd'
                pos_days = nonzeros((1:days)'.*(D(:,2) == month & wd== d_o_w & D(:,3)>14 & D(:,3)<=21));
                holiday(pos_days) = true;
            case '4th'
                pos_days = nonzeros((1:days)'.*(D(:,2) == month & wd== d_o_w & D(:,3)>21 & D(:,3)<=28));
                holiday(pos_days) = true;
        end
    end
end
end%Ends function find_holidays