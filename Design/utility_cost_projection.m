function projection = utility_cost_projection(utility,user,years,state)
% Utility is either 'gas' or 'electric'
% user is either, 'Commercial', 'Residential', or 'Industrial', if the utility is gas there is also an option for 'Electric Utility'
% years is the # of years forward to project
% state is the two letter state abreviation
global Model_dir
if strcmp(utility,'gas')
    load(fullfile(Model_dir,'System_Library', 'NatGas','RateData','GasRate'))
    Date = GasRate.Date(:,1);
    Data = GasRate.(state).(user);
    yH = floor(length(Date)/12);
    annualAvg = zeros(1,yH);
    monthAvg = zeros(1,12);
    for j = 1:1:yH
        monthAvg(1:12) = monthAvg+GasRate.(state).CityGate(12*(j-1)+1:12*j);
        annualAvg(j) = mean(GasRate.(state).CityGate(12*(j-1)+1:12*j));
    end
    monthAvg = monthAvg/yH;
elseif strcmp(utility,'electric')
    load(fullfile(Model_dir,'System_Library', 'Grid','RateData','ElecRate'))
    Date = ElecRate.Date;
    Data = ElecRate.(state).(user);
    yH = floor(length(Date)/12);
    annualAvg = zeros(1,yH);
    monthAvg = zeros(1,12);
    for j = 1:1:yH
        monthAvg(1:12) = monthAvg+ElecRate.(state).Commercial(12*(j-1)+1:12*j);
        annualAvg(j) = mean(ElecRate.(state).Commercial(12*(j-1)+1:12*j));
    end
    monthAvg = monthAvg/yH;
end

X = linspace(j+1,j+years,years);
A = polyfit(linspace(1,j,j),annualAvg,2);
Y = A(1)*X.^2+A(2)*X+A(3);
if A(2)<0 || A(1)<0
    A = polyfit(linspace(1,j,j),annualAvg,1);
    Y = A(1)*X+A(2);
    dif = Y(1)-(annualAvg(end)+A(1));
    Y = Y-dif+dif*sqrt(linspace(.001,1,years));
end

for j = 1:1:length(Y)
    projection(12*(j-1)+1:12*j) = Y(j)+monthAvg-mean(monthAvg);
end
projection = projection - (projection(1)-Data(end));