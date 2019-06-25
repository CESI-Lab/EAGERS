function forecast = arima_wsu(date,prev,options)
forecast.Timestamp = date;
freq = 1; %period of repetition (1 = 1 day)
res = options.Resolution/24;
n_o = round(freq/res)+1;
hoz = options.Horizon;
n_n = round(hoz/24/res);
future = linspace(date(1),date(1)-res + hoz,n_n)';
date2 = [hist;future];
a = 60.9/100;
b= 45.3/100;
f = fieldnames(prev);
f = f(~strcmp('Timestamp',f));
for j = 1:1:length(f) %repeat for electric, cooling, heating, and steam as necessary
    if isstruct(prev.(f{j}))
        s = fieldnames(prev.(f{j}));
        for k = 1:1:length(s)
            if isnumeric(prev.(f{j}).(s{k}))
                n_d = length(prev.(f{j}).(s{k})(1,:));
                r = zeros(n_o+n_n,n_d);
                r(1:n_o,:) = prev.(f{j}).(s{k});
                d1 = r(2:end,:) - r(1:end-1,:);%kW/timestep
                for i = n_o+1:(n_o+n_n)
                    d1(i-1,:) = a*d1(i-2,:) + b*d1(i-n_o,:);%update the delta to include the new prediction 
                    r(i,:) = d1(i-1,:) + r(i-1,:);
                end
                forecast.(f{j}).(s{k}) = interp1(date2,r,date);
            end
        end
    elseif isnumeric(prev.(f{j}))
        n_d = length(prev.(f{j})(1,:));
        r = zeros(n_o+n_n,n_d);
        r(1:n_o,:) = prev.(f{j});
        d1 = r(2:end,:) - r(1:end-1,:);%kW/timestep
        for i = n_o+1:(n_o+n_n)
            d1(i-1,:) = a*d1(i-2,:) + b*d1(i-n_o,:);%update the delta to include the new prediction 
            r(i,:) = d1(i-1,:) + r(i-1,:);
        end
        forecast.(f{j}) = interp1(date2,r,date);
    end
end
end%Ends function arima
