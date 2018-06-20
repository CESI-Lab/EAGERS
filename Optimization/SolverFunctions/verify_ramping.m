function locked = verify_ramping(gen,subnet,qp,optimal_state,dt)
%Premise, it may not be possible to follow the output from the step-by-step
%dispatch due to ramp rate limitations. The ramp rates will limit the 
%cumulative energy generation. Storage or spare capacity in other generators 
%can make up the difference without turning something on. 
%Otherwise, something must turn on.
n_g = length(gen);
n_s = length(dt);
locked = true(n_s+1,n_g);
for i = 1:1:n_g
    if ~qp.Organize.Enabled(i)
        locked(:,i)=false;
    end
end
network_names = fieldnames(qp.Organize.Balance);
%% make sure it can shut down in time from initial condition
for i = 1:1:n_g
    if qp.Organize.Dispatchable(i) ==1
        locked(optimal_state(:,i)==0,i)=false;%Identify when optimal dispatch has things off
        if optimal_state(1,i)>0 && ~all(locked(:,i))
            r = qp.Organize.Ramping(1,i)+1;
            d = optimal_state(1,i);
            t = 1;
            while d>0
                d = d - qp.b(r);
                if d>0 && ~locked(t+1,i)
                    locked(t+1,i) = true;
                end
                t = t+1;
                r = r+qp.Organize.t1ineq;
            end
        end
    end
end

%% Make sure the loss of energy due to ramping constraints at starts & stops, can be made up by other generators or storage
for net = 1:1:length(network_names)
    inc = false(n_g,1);
    stor = false(n_g,1);
    util = [];
    out = subnet.(network_names{net}).abbreviation;
    for i = 1:1:n_g
        if strcmp(gen(i).Type,'Solar') 
            %%don't include
        elseif strcmp(gen(i).Type,'Utility')
            if isfield(gen(i).QPform.output,out) || (strcmp(out,'DC') && isfield(gen(i).QPform.output,'E'))
                util = i;
            end
        elseif ismember(gen(i).Type,{'Thermal Storage';'Electric Storage';}) 
            if isfield(gen(i).QPform.output,out)
                stor(i) = true;
            end
        elseif strcmp(gen(i).Type,'Chiller') 
            if strcmp(out,'C')
                inc(i) = true;
            end
        elseif isfield(gen(i).QPform.output,out)
            inc(i) = true;
        end
    end
    if any(stor) %&& any(StoredEnergy>0)
        locked = check_capacity(gen,subnet,qp,optimal_state,stor,locked,inc,out,dt);
    else
        locked = turn_gen_on(gen,subnet,qp,optimal_state,locked,out,network_names{net},inc,util,dt);
    end
end
end%ends function verify_ramping

function locked = check_capacity(gen,subnet,qp,optimal_state,stor,locked,inc,out,dt)
%start something earlier, or leave something on, if there is not enough
%storage capacity to makeup the difference between the optimal dispatch and
%the maximum output given when the generator turns on/off
n_g = length(inc);
n_s = length(dt);
stored_energy = zeros(n_s,1);
charge_eff = zeros(n_g,1);
for i = 1:1:n_g
    if stor(i)%find the cumulative stored energy of this type
        buff = gen(i).QPform.Stor.UsableSize*(gen(i).VariableStruct.Buffer/100);
        stored_energy = stored_energy + max(0,(optimal_state(2:end,i)-buff)*gen(i).QPform.Stor.DischEff);
        charge_eff(i) = gen(i).QPform.Stor.ChargeEff;
    end
end
avg_charge_eff = mean(nonzeros(charge_eff));
[max_out,constraint] = gen_limit(gen,optimal_state,locked,dt);
if strcmp(out,'H') %overwrite electric dispatch setting with what the heat dispatch is for CHP generators
    [dispatch,heat_consumed] = heat_disp(gen,subnet,optimal_state,locked,qp,dt);
else
    dispatch = optimal_state;
end
for i = 1:1:n_g
    if inc(i)
        starts = nonzeros((1:n_s)'.*((locked(2:end,i)-locked(1:n_s,i))>0));
        k = 1;
        while ~isempty(starts) && k<=length(starts)
            t = starts(k);
            d = zeros(n_s,1);
            s = zeros(n_s,1);
            new_on = [];
            while t<=n_s && max_out.(out)(t+1,i)<dispatch(t+1,i)
                if t>starts(k)
                    d(t) = d(t-1);
                end
                d(t) = d(t)+dispatch(t+1,i)-max_out.(out)(t+1,i);%cumulative difference between desired output and feasible output
                t = t+1;
            end
            t_end = max(1,t-1);
            for t = 1:1:t_end%cumulative spare capacity in other generators
                if t>1
                    s(t) = s(t-1);
                end
                if strcmp(out,'H')
                    for j = 1:1:n_g 
                        if inc(j) && j~=i
                            s(t) = s(t) + max_out.(out)(t+1,j);
                        end
                        s(t) = s(t) - heat_consumed(t);
                    end
                else
                    for j = 1:1:n_g 
                        if inc(j) %&& j~=i
                            s(t) = s(t) + max_out.(out)(t+1,j) - dispatch(t+1,j);
                        end
                    end
                end
            end
            if any((d - s) > stored_energy) && ~isinf(constraint(starts(k),i)) && ~(constraint(starts(k),i)==0) && ~locked(constraint(starts(k),i)+1,i)
                %update Locked, MaxOut, and StoredEnergy if the generator was started earlier
                new_on = [constraint(starts(k),i), new_on];
                locked(new_on+1,i) = true;
                old_max = max_out.(out)(new_on+1,i);
                [max_out,constraint] = gen_limit(gen,optimal_state,locked,dt);
                optimal_state(new_on+1,i) = max_out.(out)(new_on+1,i);
                if strcmp(out,'H') %overwrite electric dispatch setting with what the heat dispatch is for CHP generators
                    [dispatch,heat_consumed] = heat_disp(gen,subnet,optimal_state,locked,qp,dt);
                else
                    dispatch = optimal_state;
                end
                for j = 1:1:length(new_on)
                    stored_energy(new_on(j):end) = stored_energy(new_on(j):end) + (max_out.(out)(new_on(j)+1,i)-old_max(j))*avg_charge_eff;
                end
                if starts(k)>1 && locked(starts(k),i) == false %move start and check again
                    starts(k) = starts(k)-1;
                    k = k-1;
                end
            end
            k = k+1;
        end

        stops = nonzeros((1:n_s)'.*(locked(1:n_s,i)-(locked(2:end,i))>0));
        if ~isempty(stops) && stops(1) ==1 %Initial condition ramp-down was taken care of earlier
            if length(stops)>1
                stops = stops(2:end);
            else
                stops = [];
            end
        end
        k = 1;
        while ~isempty(stops) && k<=length(stops)
            d = zeros(n_s,1);
            s = zeros(n_s,1);
            new_on = [];
            next_on = starts(starts>stops(k));
            if isempty(next_on)
                next_on = n_s;
            else
                next_on = next_on(1);
            end
            t_start_down = 1;
            while dispatch(t_start_down+1,i)<max_out.(out)(t_start_down+1,i) && t_start_down<stops(k)
                t_start_down = t_start_down+1;%find first time it exceeds ramping constraint
            end
            for t = 1:1:next_on %cumulative spare capacity in other generators
                if t>1
                    s(t) = s(t-1);
                    d(t) = d(t-1);
                end
                if t>=t_start_down
                    d(t) = max(0,d(t)+dispatch(t+1,i)-max_out.(out)(t+1,i));%difference between desired output and feasible output
                end
                
                if strcmp(out,'H')
                    for j = 1:1:n_g 
                        if inc(j) && j~=i
                            s(t) = s(t) + max_out.(out)(t+1,j);
                        end
                        s(t) = s(t) - heat_consumed(t);
                    end
                else
                    for j = 1:1:n_g 
                        if inc(j) && j~=i
                            s(t) = s(t) + max_out.(out)(t+1,j) - dispatch(t+1,j);
                        end
                    end
                end
            end
            if any((d - s) > stored_energy) && constraint(stops(k)-1,i)>0 && ~isinf(constraint(stops(k)-1,i)) && ~locked(constraint(stops(k)-1,i)+1,i) 
                new_on = [new_on,constraint(stops(k)-1,i)];
                locked(new_on+1,i) = true;
                old_max = max_out.(out)(new_on+1,i);
                [max_out,constraint] = gen_limit(gen,optimal_state,locked,dt);
                optimal_state(new_on+1,i) = max_out.(out)(new_on+1,i);
                if strcmp(out,'H') %overwrite electric dispatch setting with what the heat dispatch is for CHP generators
                    [dispatch,heat_consumed] = heat_disp(gen,subnet,optimal_state,locked,qp,dt);
                else
                    dispatch = optimal_state;
                end
                for j = 1:1:length(new_on)
                    stored_energy(new_on(j):end) = stored_energy(new_on(j):end) + (max_out.(out)(new_on(j)+1,i)-old_max(j))*avg_charge_eff;
                end
                if stops(k)<n_s && locked(stops(k)+1,i) == false %move stop and check again
                    stops(k) = stops(k)+1;
                    k = k-1;
                end
            end
            k = k+1;
        end
    end
end
end%Ends function checkStorageCapacity

function [locked,feas] = turn_gen_on(gen,subnet,qp,optimal_state,locked,out,net,inc,util,dt)
n_g = length(gen);
n_s = length(dt);
lgen = [];
net_demand = zeros(n_s,1);%% Find net demands
for n = 1:1:length(qp.Organize.Balance.(net)) %run through all the nodes in this network
    req = qp.Organize.Balance.(net)(n);%balance at this node (t = 1)
    req = req:qp.Organize.t1Balances:((n_s-1)*qp.Organize.t1Balances+req);%balance at this node (t = 1:nS)
    net_demand = net_demand + qp.beq(req);
end
for i = 1:1:n_g  %add demands from electric or absorption chillers
    if strcmp(gen(i).Type,'Chiller') && ((strcmp(out,'H') && strcmp(gen(i).Source,'Heat')) || ((strcmp(out,'E') && strcmp(gen(i).Source,'Electricity'))) )
        net_demand = net_demand + qp.constDemand.(out).load(:,i).*locked(2:end,i);
        lgen(end+1) = i;
        n = gen(i).QPform.(net).subnetNode;
        req = qp.Organize.Balance.(net)(n);%balance at this node (t = 1)
        req = req:qp.Organize.t1Balances:((n_s-1)*qp.Organize.t1Balances+req);%balance at this node (t = 1:nS)
        for t = 1:1:n_s
            if locked(t+1,i)
                C = 0;
                j = 1;
                states = qp.Organize.States{i} + (t-1)*qp.Organize.t1States;
                while j<=length(states) && C<optimal_state(t+1,i)
                    c = min(optimal_state(t+1,i) - C,qp.ub(states(j)));
                    C = C + c;
                    net_demand(t) = net_demand(t) - c*qp.Aeq(req(t),states(j));
                    j = j+1;
                end
            end
        end
    end
end
%%No storage to check 
% for i = 1:1:nG  %add demands for charging/discharging storage
%     if any(strcmp(Gen(i).Type, {'Thermal Storage';'Electric Storage'})) && isfield(Gen(i).QPform.output,out)
%         for t = 1:1:nS
%             if OptimalState(t,i)<OptimalState(t+1,i)%charging
%                 NetDemand(t) = NetDemand(t) + (OptimalState(t+1,i) - OptimalState(t,i))*dt(t)/Gen(i).QPform.Stor.ChargeEff;
%             else
%                 NetDemand(t) = NetDemand(t) - (OptimalState(t,i) - OptimalState(t+1,i))*dt(t)*Gen(i).QPform.Stor.DischEff;
%             end
%         end 
%     end
% end
feas = false;
while ~feas
    [max_out,constraint] = gen_limit(gen,optimal_state,locked,dt);
    if strcmp(out,'H') %overwrite electric dispatch setting with what the heat dispatch is for CHP generators
        [~,net_demand] = heat_disp(gen,subnet,optimal_state,locked,qp,dt);
    end
    constrained_gen = sum(max_out.(out)(2:end,inc),2);
    if ~isempty(util)
        constrained_gen = constrained_gen + max_out.(out)(2:end,util);
    end
    if all(constrained_gen>=(net_demand-1e-9))%fix due to rounding errors
        feas = true;
    else
        feas = false;
        t = 1;
        while sum(max_out.(out)(t+1,inc))>=(net_demand(t)-1e-9)%fix due to rounding errors
            t = t+1;
        end
        diff = zeros(n_g,1);
        if any(constraint(t,inc)>0 & constraint(t,inc)<n_s+1) %at leat one generator in this category is constrained by a startup or shutdown at this time
            %select apropriate one to start early or keep on later
            for i = 1:1:n_g
                if inc(i)
                    if constraint(t,i)>0 && constraint(t,i)<n_s+1
                        diff(i) = abs(constraint(t,i)-t);
                    else
                        diff(i) = n_s+1;
                    end
                else
                    diff(i) = inf;
                end
            end
            [~,I] = min(diff);
            locked(constraint(t,I)+1,I) = true;
        elseif any(~locked(t+1,inc))
            %turn on something else that has the same type of output
            for i = 1:1:n_g
                if inc(i)
                    if ~locked(t+1,i)
                        if t==1
                            diff(i) = n_s-t - max((n_s-t:-1:1)'.*locked(t+2:n_s+1,i));
                        elseif t == n_s
                            diff(i) = t - max((1:t-1)'.*locked(2:t,i));
                        else
                            diff(i) = min(t - max((1:t-1)'.*locked(2:t,i)),n_s-t - max((n_s-t:-1:1)'.*locked(t+2:n_s+1,i)));
                        end
                    else
                        diff(i) = n_s+1;
                    end
                else
                    diff(i) = inf;
                end
            end
            [~,I] = min(diff);
            locked(t+1,I) = true;
        elseif ~isempty(lgen) && any(locked(t+1,lgen))
            %turn off something that is a load (e.g. a absorption chiller when balancing the heat
            for i = 1:1:length(lgen)
                if locked(t+1,lgen(i))
                    locked(t+1,lgen(i)) = false;
                end
            end
        else
            feas = true;%nothing else it can change by this logic
%             disp('Error in CheckRampRates: cant find a generator to turn on to make feasible')
        end
    end
end
end%Ends function turn_something_on

function [dispatch,heat_consumed] = heat_disp(gen,subnet,optimal_state,locked,qp,dt)
n_g = length(gen);
[n,~] = size(optimal_state);
n_s = n - 1;
dispatch = optimal_state; 
heat_consumed = zeros(n_s,1);
for n = 1:1:length(subnet.DistrictHeat.nodes) %run through all the nodes in this network
    req = qp.Organize.Balance.DistrictHeat(n);%balance at this node (t = 1)
    req = req:qp.Organize.t1Balances:((n_s-1)*qp.Organize.t1Balances+req);%balance at this node (t = 1:nS)
    load = nonzeros(qp.Organize.Demand.DistrictHeat{n}); %loads at this node
    if ~isempty(load) %need this in case there is no field Forecast.Demand
        heat_consumed = heat_consumed + qp.beq(req); %multiple demands can be at the same node, or none
    end
end
for i = 1:1:n_g
    if strcmp(gen(i).Type, 'CHP Generator')
        for t = 1:1:n_s
            cum_out = 0;
            j = 1;
            if locked(t+1,i)
                dispatch(t+1,i) = -gen(i).QPform.constDemand.H;               
            else
                dispatch(t+1,i) = 0;
            end
            states = qp.Organize.States{i} + (t-1)*qp.Organize.t1States;
            while j<=length(states)
                d = min(optimal_state(t+1,i) - cum_out,qp.ub(states(j)));
                cum_out = cum_out + d;
                dispatch(t+1,i) = dispatch(t+1,i) + d*qp.Aeq(req(t),states(j));
                j = j+1;
            end
        end
    elseif strcmp(gen(i).Type, 'Heater')

    elseif strcmp(gen(i).Type, 'Chiller') && strcmp(gen(i).Source,'Heat')
        dispatch(:,i) = 0;
        heat_consumed = heat_consumed + qp.constDemand.H.load(:,i).*locked(2:end,i);
        for t = 1:1:n_s
            if locked(t+1,i)
                cum_out = 0;
                j = 1;
                states = qp.Organize.States{i} + (t-1)*qp.Organize.t1States;
                while j<=length(states) && cum_out<optimal_state(t+1,i)
                    d = min(optimal_state(t+1,i) - cum_out,qp.ub(states(j)));
                    cum_out = cum_out + d;
                    heat_consumed(t) = heat_consumed(t) - d*qp.Aeq(req(t),states(j));
                    j = j+1;
                end
            end
        end
    elseif strcmp(gen(i).Type, 'Thermal Storage')
        if strcmp(gen(i).Source,'Heat')
            for t = 1:1:n_s
                if optimal_state(t,i)<optimal_state(t+1,i)%charging
                    heat_consumed(t) = heat_consumed(t) + (optimal_state(t+1,i) - optimal_state(t,i))*dt(t)/gen(i).QPform.Stor.ChargeEff;
                else
                    heat_consumed(t) = heat_consumed(t) - (optimal_state(t,i) - optimal_state(t+1,i))*dt(t)*gen(i).QPform.Stor.DischEff;
                end
            end  
        else
            dispatch(:,i) = 0;
        end
    else
        dispatch(:,i) = 0;
    end
end
end%Ends function heat_disp