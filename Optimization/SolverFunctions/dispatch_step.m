function [gen_output,flag] = dispatch_step(gen,building,fluid_loop,subnet,options,one_step,date,forecast,scale_cost,dt,first_profile)
%DISPATCH_STEP
%
% stor_power is the amount of power coming from (positive) or going into (negative) the storage device at each timestep according to the first dispatch
%
% Flag values:
%   0 -- Standard operation.
%   1 -- No feasible combinations at one or more time steps.

% Initialize flag.
flag = 0;

n_g = length(gen);
n_b = length(building);
n_l = length(one_step.Organize.Transmission);
n_fl = length(fluid_loop);
temperatures.build = zeros(2,n_b);
for i = 1:1:n_b
    temperatures.build(1,i) = building(i).Tzone;
    temperatures.build(2,i) = building(i).Twall;
end
temperatures.fluid_loop = zeros(1,n_fl);
for i = 1:1:n_fl
    temperatures.fluid_loop(i) = fluid_loop(i).fluid_temperature;
end

limit = 'initially constrained';
start_cost = zeros(1,n_g);
ic = zeros(1,n_g);
for i = 1:1:n_g
    if isfield(gen(i).VariableStruct, 'StartCost')
        start_cost(i) = gen(i).VariableStruct.StartCost;
    end
    ic(i) = gen(i).CurrentState(1);
end
[n_s,~] = size(scale_cost);
gen_output = zeros(n_s+1,n_g+n_l+n_b+n_fl);
gen_output(1,1:n_g) = ic;


stor_power = zeros(n_s,n_g);
i_best = ones(n_s,1);
cost = cell(n_s,1);
verified = cell(n_s,1);
qp = cell(n_s,1);
binary_comb = cell(n_s,1);
disp_comb = cell(n_s,1);
for t = 1:1:n_s %for every timestep
    net_demand = agregate_demand(forecast,t);
    dem_h = 0;
    if (isfield(net_demand,'H'))
        dem_h = net_demand.H;
    end
    v_h = value_heat(gen,[ic;first_profile(t+1,:)],dem_h,dt(t));
    stor_power(t,:) = find_stor_power(gen,ic,first_profile(t+1,:),dt(t));
    for i = 1:1:n_g
        if strcmp(gen(i).Type,'Hydro Storage')
            stor_power(t,i) = first_profile(t+1,i) + (first_profile(t,i) - ic(i));
        end
    end
    marginal = update_mc(gen,first_profile(t+1,:),scale_cost(t,:),dt(t),v_h);%update marginal cost
    qp{t} = update_matrices_step(gen,building,fluid_loop,subnet,options,one_step,forecast,scale_cost(t,:),marginal,stor_power(t,:),dt,ic,first_profile,limit,t,temperatures);
    
    [verified{t},cost{t},binary_comb{t},disp_comb{t}] = test_min_cases(qp{t},gen,options,net_demand,marginal,scale_cost(t,:),first_profile(t+1,:),dt(t));   
    if isempty(verified{t})
        flag = 1;
        gen_output = first_profile;
%         disp(strcat('no feasible combinations to test at step', num2str(t)))
        return
    end
    %update Initial conditions and building temperatures
    [~,i_best(t)] = min(cost{t});
    best_dispatch = disp_comb{t}(i_best(t),:);
    [gen_output(t+1,:),ic] = update_ic(gen,ic,best_dispatch,first_profile(t+1,:),dt(t),limit);%only updates the storage states in IC
    temperatures = buildings_step(building,forecast,temperatures,qp{t},binary_comb{t}(i_best(t),:),dt,t);
    for i = 1:1:n_fl
        temperatures.fluid_loop(i) = fluid_loop_step(gen,subnet.CoolingWater.Equipment{i},fluid_loop(i),temperatures.fluid_loop(i),best_dispatch,dt(t)*3600);
    end
end
% %change best dispatch based on re-start costs
gen_output = minimize_start_costs(qp,gen,[date;forecast.Timestamp],gen_output,stor_power,binary_comb,disp_comb,cost,verified,start_cost,dt);
if isfield(forecast,'Renewable')
    gen_output(2:n_s+1,1:n_g) = gen_output(2:n_s+1,1:n_g) + forecast.Renewable;
end
end% Ends function dispatch_step

function [gen_output,ic] = update_ic(gen,ic,best_dispatch,first_profile,dt,limit)
n_g = length(gen);
gen_output = best_dispatch;
if strcmp(limit,'unconstrained')
    ic = [];
else
    ec = best_dispatch(1:n_g);
    if ~isempty(first_profile)
        for i = 1:1:n_g
            if ismember(gen(i).Type,{'Electric Storage';'Thermal Storage'})
                ec(i) = first_profile(i) - best_dispatch(i)*dt/gen(i).QPform.Stor.DischEff;
                if ec(i)<0
                    disp(strcat('Warning: ',gen(i).Name,'_ is going negative by_',num2str(ec(i)/gen(i).QPform.Stor.UsableSize*100),'%'))
                    ec(i) = 0;
                elseif ec(i)>gen(i).QPform.Stor.UsableSize
                    disp(strcat('Warning: ',gen(i).Name,'_ is exceeding max charge by_',num2str((ec(i)-gen(i).QPform.Stor.UsableSize)/gen(i).QPform.Stor.UsableSize*100),'%'))
                    ec(i) = gen(i).QPform.Stor.UsableSize;                    
                end
            end
        end
    end
    gen_output(1,1:n_g) = ec;
    if strcmp(limit,'constrained')%if its constrained but not initially constrained then make the last output the initial condition
        ic = ec;
    else %if strcmp(limit,'initially constrained')
        for i = 1:1:n_g
            switch gen(i).Type
                case {'Electric Storage';'Thermal Storage';'Hydrogen Storage';}
                    ic(i) = ec(i); %update the state of storage
                case 'Hydro Storage'
                    %% update something. Need better calculation of min/max power at next step according to remaining state of reservoir
                    ic(i) = gen_output(1,i);
            end
        end
    end
end
end%ends function update_ic

function net_demand = agregate_demand(forecast,t)
net_demand = [];
if isfield(forecast,'Demand')
    outs = fieldnames(forecast.Demand);
    for j = 1:1:length(outs)
        net_demand.(outs{j}) = sum(forecast.Demand.(outs{j})(t,:));
    end
end
if isfield(forecast,'Building')
    outs2 = {'E';'H';'C';};
    for j = 1:1:length(outs2)
        if ~isfield(net_demand,outs2{j})
            net_demand.(outs2{j}) = 0;
        end
        net_demand.(outs2{j}) = net_demand.(outs2{j}) + sum(forecast.Building.(strcat(outs2{j},'0'))(t,:));
    end
end
end%Ends function agregate_demand

function v_h = value_heat(gen,dispatch,excess_heat,dt)
n_g = length(gen);
for i = 1:1:n_g
    if strcmp(gen(i).Type,'Thermal Storage') && strcmp(gen(i).Source,'Heat')
        loss = dt*(gen(i).QPform.Stor.SelfDischarge*gen(i).QPform.Stor.UsableSize);
        d_SOC = dispatch(2,i) - dispatch(1,i);%expected output of storage in kWh to reach the SOC from the 1st dispatch (penalties are always pushing it back on track if it needed more storage than expected somewhere)
        if (d_SOC/dt + loss)<0 %discharging
            excess_heat = excess_heat -(d_SOC/dt + loss)*gen(i).QPform.Stor.DischEff;
        else %charging
            excess_heat = excess_heat -(d_SOC/dt +loss)/gen(i).QPform.Stor.ChargeEff; 
        end
    end
    if strcmp(gen(i).Type,'CHP Generator')
        j = 1;
        d = 0;
        states = gen(i).QPform.states(1:nnz(~cellfun('isempty',gen(i).QPform.states(:,end))),end);
        if dispatch(2,i)>0
            excess_heat = excess_heat - gen(i).QPform.constDemand.H;%for CHP generators, this constDemand is negative
        end
        while j<=length(states) && dispatch(2,i)-d>0
            x_i = min(dispatch(2,i)-d,gen(i).QPform.(states{j}).ub(end));
            excess_heat = excess_heat + x_i*gen(i).QPform.output.H(min(j,length(gen(i).QPform.output.H(:,1))),1);
            d = d + x_i;
            j = j+1;
        end
    end
end
v_h = excess_heat<=1e-1;
end%Ends function value_heat

function stor_power = find_stor_power(gen,ic,dispatch,dt)
n_g = length(gen);
stor_power = zeros(1,n_g);  
for i = 1:1:n_g
    if any(strcmp(gen(i).Type,{'Electric Storage';'Thermal Storage';'Hydrogen Storage';}))
        loss = (gen(i).QPform.Stor.SelfDischarge*gen(i).QPform.Stor.UsableSize);
        d_SOC = dispatch(1,i) - ic(i);%expected output of storage in kWh to reach the SOC from the 1st dispatch (penalties are always pushing it back on track if it needed more storage than expected somewhere)
        if (d_SOC/dt + loss)<0 %discharging
            stor_power(1,i) = -(d_SOC/dt + loss)*gen(i).QPform.Stor.DischEff;
        else %charging
            stor_power(1,i) = -(d_SOC/dt +loss)/gen(i).QPform.Stor.ChargeEff; 
        end
    end
end
end%Ends function find_stor_power


function temperatures = buildings_step(building,forecast,temperatures,qp,binary_comb,dt,t)
n_b = length(building);
if n_b>0%%re-do best cost case
    qp_test = disable_generators_step(qp,binary_comb);%Disable generators here
    x = call_solver(qp_test);
    [~,~,~,~,~,~,heating,cooling] = sort_solution_step(x,qp_test);
%             [best_dispatch,line_loss(t,:),excess_heat(t,:),excess_cool(t,:),hydro_soc(t,:),temperature(t,:),heating(t,:),cooling(t,:)] = sort_solution_step(x,qp_test);
    for i = 1:1:n_b
        [zone,wall] = building_simulate(building(i),...
            forecast.Weather.Tdb(t),forecast.Weather.RH(t),dt(t)*3600,...
            forecast.Building.InternalGains(t,i),forecast.Building.ExternalGains(t,i),...
            cooling(1,i),heating(1,i),forecast.Building.AirFlow(t,i),...
            forecast.Building.Damper(t,i),temperatures.build(1,i),temperatures.build(2,i));
        temperatures.build(:,i) = [zone(2);wall(2)];
    end
end
end%Ends function buildings_step
