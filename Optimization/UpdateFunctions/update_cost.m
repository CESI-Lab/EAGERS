function scale_cost = update_cost(timestamp,gen)
n_s = length(timestamp);
n_g = length(gen);
source = cell(n_g,1);
utility_rate = zeros(n_s,n_g);
for i = 1:1:n_g
    if strcmp(gen(i).Type,'Utility')
        source(i) = {gen(i).Source};
        utility = gen(i).VariableStruct;
        if strcmp(source(i), 'Electricity')
            day = weekday(timestamp);
            [y,m,d] = datevec(timestamp(1));
            win_start = datenum([y,utility.WinStartMonth,utility.WinStartDay]);
            if m>utility.WinStartMonth || (m==utility.WinStartMonth && d>=utility.WinStartDay)
                sum_start = datenum([y+1,utility.SumStartMonth,utility.SumStartDay]);
            else
                sum_start = datenum([y,utility.SumStartMonth,utility.SumStartDay]);
            end
            rate = zeros(n_s,1);
            [~,~,~,h,~,~] = datevec(timestamp);
            h(h==0) = 24; %zeroth hour is column 1 in rate matrix
            for t = 1:1:length(h) 
                if timestamp(t)>=sum_start && timestamp(t)<win_start
                    rate(t) = utility.SumRates(utility.SumRateTable(day(t),h(t),1));
                else
                    rate(t) = utility.WinRates(utility.WinRateTable(day(t),h(t),1));
                end
            end
            utility_rate(:,i) = rate;%rate in $/kWh
        elseif strcmp(source(i), 'NG')|| strcmp(source(i), 'Diesel')
            %gas utility set up as daily prices for a leap year
            date1 = datevec(utility.Timestamp(1));
            utility.Timestamp(end+1) = utility.Timestamp(end)+1;
            utility.Rate(end+1) = utility.Rate(1);
            year1 = date1(1);
            datenow = datevec(timestamp);
            interp_date = datenum([year1*ones(length(timestamp),1), datenow(:,2:end)]);
            utility_rate(:,i) = interp1(utility.Timestamp,utility.Rate,interp_date)/293.1; %interpolate & convert gas rate from $/MMBTU to $/kWh;
        else
            %% need to add something for district heating/cooling
        end
    end
end

scale_cost = zeros(n_s,n_g);
for i = 1:1:length(gen)
    if strcmp(gen(i).Type,'Utility')
        scale_cost(:,i) = utility_rate(:,i);
    elseif isempty(strfind(gen(i).Type,'Storage')) && ~strcmp(gen(i).Source, 'Renewable') && ~strcmp(gen(i).Type, 'Hydro') %not storage or renewable or hydro
        util_index = nonzeros((1:n_g)'.*strcmp(gen(i).Source,source));
        if isempty(util_index)
            scale_cost(:,i) = 1; %no utility (don't scale costs)
        else
            scale_cost(:,i) = utility_rate(:,util_index);
        end
    end
end