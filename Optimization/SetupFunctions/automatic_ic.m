function [gen,fluid_loop,flag] = automatic_ic(gen,building,fluid_loop,subnet,date,one_step,options,data_t0)

% Flag values:
%   0 -- Standard operation.
%   1 -- No feasible combination.

%determine the net energy demands in each category
n_g = length(gen);
net_demand = [];
if isfield(data_t0,'Demand')
    outs = fieldnames(data_t0.Demand);
    for j = 1:1:length(outs)
        net_demand.(outs{j}) = sum(data_t0.Demand.(outs{j}),2);
    end
end

temperatures.build = [];
n_b = length(building);
for i = 1:1:n_b
    temperatures.build(i,:) = data_t0.Building(i).Tzone;
    Outs2 = {'E';'H';'C';};
    for j = 1:1:length(Outs2)
        if ~isfield(net_demand,Outs2{j})
            net_demand.(Outs2{j}) = sum(data_t0.Building.(strcat(Outs2{j},'0')));
        else
            net_demand.(Outs2{j}) = net_demand.(Outs2{j}) + sum(data_t0.Building.(strcat(Outs2{j},'0')));
        end
    end
end

n_fl = length(fluid_loop);
for i = 1:1:n_fl
    fluid_loop(i).fluid_temperature = fluid_loop(i).nominal_return_temperature;
    temperatures.fluid_loop(i) = fluid_loop(i).fluid_temperature;
    %set upper and lower bound equal to this temperature to initialize at steady state
    n_l = length(one_step.Organize.Transmission);
    state = one_step.Organize.States{n_g+n_l+n_b+i};
    one_step.lb(state) = fluid_loop(i).fluid_temperature-.01;
    one_step.ub(state) = fluid_loop(i).fluid_temperature;
end

scale_cost = update_cost(date,gen);%% All costs were assumed to be 1 when building matrices, update Generator costs for the given time

if isfield(data_t0,'Weather') && isfield(data_t0.Weather,'DNIWm2')
    data_t0.Renewable = renewable_output(gen,subnet,date,data_t0.Weather.DNIWm2);
end

marginal = update_mc(gen,[],scale_cost,[],0);%update marginal cost
QP = update_matrices_step(gen,building,fluid_loop,subnet,options,one_step,data_t0,scale_cost,marginal,[],options.Resolution,[],[],'unconstrained',1,temperatures);

[~,cost,~,disp_comb] = test_min_cases(QP,gen,options,net_demand,marginal,scale_cost,[],[]);   
if isempty(disp_comb)
    flag = 1;
else
    flag = 0;
    [~,i_best] = min(cost);
    ic = disp_comb(i_best,:);
    gen = set_ic(gen,ic,50);
end
end%Ends function automatic_ic

function gen = set_ic(gen,ic,stor_perc)
n_g = length(gen);
lb = zeros(1,n_g);
for i=1:1:n_g
    switch gen(i).Type
        case {'Electric Generator';'CHP Generator';'Chiller';'Heater';'Electolyzer';'Hydrogen Generator';'Cooling Tower'}
            states = gen(i).QPform.states(1:nnz(~cellfun('isempty',gen(i).QPform.states(:,end))),end);
            for j = 1:1:length(states)
                lb(i) = lb(i) + gen(i).QPform.(states{j}).lb(2);
            end
        case {'Electric Storage';'Thermal Storage';'Hydrogen Storage';}
            ic(i) = stor_perc/100*gen(i).QPform.Stor.UsableSize; % IC = halfway charged energy storage
    end
    gen(i).CurrentState(1) = ic(i);
    gen(i).Status = ic(i)>lb(i);
end
end%ends function set_ic