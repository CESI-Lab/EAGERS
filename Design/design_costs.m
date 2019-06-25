function [cost, npc, scaled_monthly_costs] = design_costs(gen, design, options, years, equip_costs, discount_rate, dispatchFlag)
%DESIGN_COSTS

% Set infinite cost for infeasible dispatch.
% dispatchFlag:
%   0 -- feasible
%   1, 2, 3 -- infeasible
if dispatchFlag > 0
    cost = Inf;
    npc = Inf;
    scaled_monthly_costs = Inf;
    return
end

% Rest of function.
gen_disp = design.GeneratorState;
fuel_input = zeros(size(gen_disp));

n_g = length(gen);
electric_utility = [];
for i = 1:1:n_g
    skip = false;
    if strcmp(gen(i).Type,'Electric Generator') || strcmp(gen(i).Type,'CHP Generator') || strcmp(gen(i).Type,'Hydrogen Generator')
        if isfield(gen(i).Output,'Electricity')
            eff = gen(i).Output.Electricity;
        elseif isfield(gen(i).Output,'DirectCurrent')
            eff = gen(i).Output.DirectCurrent;
        end
        cap = gen(i).Output.Capacity*gen(i).Size;
    elseif strcmp(gen(i).Type,'Electrolyzer')
        eff = gen(i).Output.Hydrogen;
        cap = gen(i).Output.Capacity*gen(i).Size;
    elseif strcmp(gen(i).Type,'Heater')
        eff = gen(i).Output.Heat;
        cap = gen(i).Output.Capacity*gen(i).Size;
    elseif strcmp(gen(i).Type,'Utility') && strcmp(gen(i).Source,'Electricity')
        electric_use = gen_disp(:,i);
        electric_utility = i;
        skip = true;
    else
        skip = true;%skip chillers whose input shows up as electric or heat generation, skip storage devices
    end
    if ~skip %dont add the cost of a chiller if you ran E and C simultaneously, or you will double count the chiller demand
        gen_eff = interp1(cap,eff,gen_disp(:,i));
        fuel_input(:,i) = gen_disp(:,i)./gen_eff;
        fuel_input(gen_disp(:,i)==0,i) = 0;
    end
end

scale_cost = update_cost(design.Timestamp,gen)*options.Resolution;%% All costs were assumed to be 1 when building matrices
fuel_cost = sum(fuel_input.*scale_cost,2);
if ~isempty(electric_utility) && gen(electric_utility).VariableStruct.MinImportThresh<0
    if gen(electric_utility).VariableStruct.SellBackRate>0
        sell_back = gen(electric_utility).VariableStruct.SellBackRate*options.Resolution;
    else
        sell_back = gen(electric_utility).VariableStruct.SellBackPerc*scale_cost(:,electric_utility);
    end
else
    sell_back = 0; %no sellback allowed
end
elec_use_cost = max(0,electric_use).*scale_cost(:,electric_utility);
elec_sell_back = -min(0,electric_use).*sell_back;

n = length(equip_costs);
operations_and_maintenance = zeros(n,1);
finance_monthly = zeros(n,1);
capital_cost = zeros(n,1);
for i = 1:1:n
    if ~isempty(equip_costs(i).Cost)
        capital_cost(i) = equip_costs(i).Cost;
        operations_and_maintenance(i) = equip_costs(i).OandM; %yearly maintenance cost
        loan_amount = equip_costs(i).Cost*equip_costs(i).Financed/100; %principle loan amount
        monthly_rate = equip_costs(i).LoanRate/12/100; 
        loan_term = equip_costs(i).LoanTerm;
        finance_monthly(i)=(monthly_rate*loan_amount/(1-(1+monthly_rate)^(-loan_term*12))); % c*L/(1-(1+c)^(-n)) is the same as L*(c*(1+c)^n)/((1+c)^n -1)
    end
end

monthly_costs = zeros(0,6); 
d = datevec(design.Timestamp(1));
m = 0;
frac_of_month = (datenum([d(1) d(2)+1 1]) - design.Timestamp(1))/(datenum([d(1) d(2)+1 1])-datenum([d(1) d(2) 1]));
while datenum([d(1) d(2)+m 1])<design.Timestamp(end)
    index = max(1,nnz(design.Timestamp<=datenum([d(1) d(2)+m 1]))):1:nnz(design.Timestamp<datenum([d(1) d(2)+m+1 1]));
    m = m+1;
    startup_cost = zeros(n_g,1);
    for i = 1:1:n_g
        if isfield(gen(i).VariableStruct, 'StartCost')
            nStartups = nnz((gen_disp(index(2:end),i)>0) & (gen_disp(index(1:end-1),i)==0));
            startup_cost(i) = nStartups*gen(i).VariableStruct.StartCost;
        end
    end
    monthly_costs(end+1,1) = sum(finance_monthly);
    monthly_costs(m,2) = sum(operations_and_maintenance)/12;
    monthly_costs(m,3) = sum(startup_cost);
    monthly_costs(m,4) = demand_charge(electric_use(index),design.Timestamp(index),gen(electric_utility).VariableStruct);
    monthly_costs(m,5) = sum(elec_use_cost(index))/frac_of_month - sum(elec_sell_back(index))/frac_of_month;
    monthly_costs(m,6) = sum(fuel_cost(index))/frac_of_month;
    frac_of_month = min(1,(design.Timestamp(end) - datenum([d(1) d(2)+m 1]))/(datenum([d(1) d(2)+m 1])-datenum([d(1) d(2)+m-1 1])));
end
cost = monthly_costs;
%find net present cost
months = 12*years;
escalator = ones(months,6);
escalator_electric = utility_cost_projection('electric','Commercial',years,'WA')';
escalator_electric = escalator_electric/escalator_electric(1);
escalator_gas = utility_cost_projection('gas','Commercial',years,'WA')';
escalator_gas = escalator_gas/escalator_gas(1);
escalator(:,4) = escalator_electric;
escalator(:,5) = escalator_electric;
escalator(:,6) = escalator_gas;
scaled_monthly_costs = zeros(months,1);
scaled_monthly_costs(1:m) = sum(monthly_costs(1:m,2:end),2);
scaled_monthly_costs(m+1:12) = sum(scaled_monthly_costs(1:m))/m;
for n = m+1:1:12%make sure there are at least 12 months
    monthly_costs(end+1,:) = sum(monthly_costs(1:m,:),1)/m;
end
for i = m+1:1:months
    scaled_monthly_costs(i) = sum(monthly_costs(rem(i-1,12)+1,2:end).*escalator(i,2:end));%monthly costs without financing charges
end
npc = net_pres_cost((1+discount_rate/100),scaled_monthly_costs+sum(finance_monthly));
scaled_monthly_costs(1) = scaled_monthly_costs(1)+sum(capital_cost);%put all capital cost in 1st month


function cost = demand_charge(use,timestamp,utility)
%find the highest demand charge of the month
day = weekday(timestamp);
[Y,~,~,~,~,~] = datevec(timestamp(1));
sum_start = datenum([Y(1),utility.SumStartMonth,utility.SumStartDay]);
win_start = datenum([Y(1),utility.WinStartMonth,utility.WinStartDay]);
if timestamp(1)>sum_start
    sum_start = datenum([Y(1)+1,utility.SumStartMonth,utility.SumStartDay]); %add a year
end
if timestamp(1)>win_start
    win_start = datenum([Y(1)+1,utility.WinStartMonth,utility.WinStartDay]);% add a year
end
if win_start<sum_start
    season = 'Sum';
else
    season = 'Win';
end
rates = utility.(strcat(season,'Rates'))(:,2);
[~,~,~,H,~,~] = datevec(timestamp);
on_peak = zeros(length(day),1);
mid_peak = zeros(length(day),1);
off_peak = zeros(length(day),1);
for t = 1:1:length(H) 
    j = utility.(strcat(season,'RateTable'))(day(t),H(t)+1);
    if j==1
        on_peak(t) = use(t);
    elseif j==2
        mid_peak(t) = use(t);
    else
        off_peak(t) = use(t);
    end
end
cost = max(on_peak)*rates(1);
cost = cost + max(mid_peak)*rates(2);
cost = cost + max(off_peak)*rates(3);

function value= net_pres_cost(irr,cashflows)
value=0;
for i=1:length(cashflows)
    value=value+cashflows(i)/(irr)^(i/12);
end
