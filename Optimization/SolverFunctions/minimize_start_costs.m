function [gen_output,i_best] = minimize_start_costs(qp,gen,date,gen_output,stor_power,binary_comb,disp_comb,cost_comb,verified,start_cost,dt)
%This function finds the n shortest segments of time the generator is on or
%off, and compares the start-up cost to the alternative configuration that
%avoids the start or re-start
% GenOutput is the bet dispatch at each time without considering start-up costs
% Alt is all of the other feasible combinations tested
% Binary is the current best on/off configuration
% StartCost is the startup cost of each generator
% dt is the duration of each time segment
% Type specifies the category of generator (currently this function is only interested in electric and CHP generators.
% include is what type of generators to consider (currently this is always electric and CHP generators)
% n is # of segments to check
n_s = length(dt);
n_g = length(gen);
include = {'Electric Generator';'CHP Generator';'Heater';'Electrolyzer';'Hydrogen Generator';'Cooling Tower';'Chiller'};
skip_on = [];
skip_off = [];
locked = true(n_s+1,n_g);
locked(1,:) = gen_output(1,1:n_g)~=0;
for t = 2:1:n_s+1
    locked(t,qp{t-1}.Organize.Dispatchable) = gen_output(t,qp{t-1}.Organize.Dispatchable)~=0;
end
inc = false(1,n_g);
for i = 1:1:n_g
    if ismember(gen(i).Type,include)
        inc(i) = true;
    end
end
i_best = ones(n_s,1);
for t = 1:1:n_s
    [~,i_best(t)] = min(cost_comb{t});
end
%find i_best(t), index that results in lowest cost once start-up is considered.
seg = 1;
on_seg = 1;
while ~isempty(on_seg) && seg<40
    %%Try and remove the shortest generator on segment, and if equal lengths, highest start cost, if same generator then first segment
    [on_seg, off_seg] = seg_length(locked,start_cost,dt,skip_on,skip_off,inc);
    if isempty(on_seg) && isempty(off_seg)
        break %no  segments to check
    elseif ~isempty(on_seg)
        [~,i_on] = min(on_seg(:,4)-start_cost(on_seg(:,1))'/(2*max(start_cost))+on_seg(:,2)/n_s); %shortest segment, and if equal lengths, highest start cost, if same generator then first segment
        k = on_seg(i_on,1); %generator index
        t1 = on_seg(i_on,2); %index when generator turns on
        t2 = on_seg(i_on,3);%index when generator shuts off
        %% Find the cheapest feasible alternative dispatch without this generator (only use generators that were on previously or will be on)
        %First try alternate generators at same time step that may have cheaper startup (or that are already on)
        [i_best,locked,no_alt,disp_comb,cost_comb,verified] = alt_generation(qp,gen,i_best,k,t1,t2,start_cost,binary_comb,disp_comb,cost_comb,verified,gen_output,locked,dt);
        if no_alt
            %Second, can it get close to the final storage capacity without this on-segment?
            [binary_comb,disp_comb,cost_comb,verified,i_best,locked,must_replace] = avoid_generation(qp,gen,date,i_best,k,t1,t2,start_cost,binary_comb,disp_comb,cost_comb,verified,gen_output,locked,dt);
            if must_replace
                %Third try moving generation earlier or later with same generator (if there is storage)
                [binary_comb,disp_comb,cost_comb,verified,i_best,locked,cant_move] = move_generation(qp,gen,i_best,k,t1,t2,binary_comb,disp_comb,cost_comb,verified,gen_output,locked,inc,dt);
                if cant_move
                    skip_on(end+1,1:4) = on_seg(i_on,:); %add to list of segments to avoid
                end
            end
        end
        gen_output = update_storage(gen,gen_output,disp_comb,i_best,stor_power,dt);
    end
    %%Try and remove the shortest generator off segment, if it is a shorter segment than the next shortest on segment
    [on_seg, off_seg] = seg_length(locked,start_cost,dt,skip_on,skip_off,inc);
    if isempty(on_seg) && isempty(off_seg)
        break %no  segments to check
    elseif ~isempty(off_seg) && (isempty(on_seg) || (min(off_seg(:,4)-start_cost(off_seg(:,1))'/(2*max(start_cost))+off_seg(:,2)/n_s) < min(on_seg(:,4)-start_cost(on_seg(:,1))'/(2*max(start_cost))+on_seg(:,2)/n_s)))%preference is to turn things off, rather than turn things on
        %third leave generator on, maybe turn off a different generator
        [~,Ioff] = min(off_seg(:,4)-start_cost(off_seg(:,1))'/(2*max(start_cost))+off_seg(:,2)/n_s); %shortest segment, and if equal lengths, highest start cost
        k = off_seg(Ioff,1); %generator index
        t1 = off_seg(Ioff,2); %index when generator turns on
        t2 = off_seg(Ioff,3);%index when generator shuts off
        [i_best,locked,cant_keep_on,disp_comb,cost_comb,verified] = leave_gen_on(qp,gen,i_best,k,t1,t2,start_cost,binary_comb,disp_comb,cost_comb,verified,gen_output,locked,dt);
        if cant_keep_on
            skip_off(end+1,1:4) = off_seg(Ioff,:); %add to list of segments to avoid
        end
        gen_output = update_storage(gen,gen_output,disp_comb,i_best,stor_power,dt);
    end 
    seg = seg+1;
end
end% Ends function minimize_start_costs

function gen_output = update_storage(gen,gen_output,disp_comb,i_best,stor_power,dt)
%pull the corresponding best dispatches with the start-cost considered
status = gen_output(1,:);
n_g = length(gen);
n_s = length(dt);
for t = 1:1:n_s
    new_status = disp_comb{t}(i_best(t),:);
    for i = 1:1:n_g
        if ismember(gen(i).Type,{'Electric Storage';'Thermal Storage'})
            loss = (gen(i).QPform.Stor.SelfDischarge*gen(i).QPform.Stor.UsableSize*gen(i).QPform.Stor.DischEff);
            if stor_power(t,i)>0 %discharging
                d_soc = (-stor_power(t,i)/gen(i).QPform.Stor.DischEff - loss)*dt(t);
            else %charging
                d_soc = (-stor_power(t,i)*gen(i).QPform.Stor.ChargeEff - loss)*dt(t);
            end
            new_d_soc = -new_status(i)*dt(t)/gen(i).QPform.Stor.DischEff;
            new_status(i) = status(i) + d_soc + new_d_soc;
            if new_status(i)>gen(i).QPform.Stor.UsableSize
                %Changed locked, but there was slack in other timesteps so storage will not actually overcharge
%                     disp(strcat('Warning, Potentially Over Charging Storage in ',num2str(sum(dt(1:t))),'_hours'))
                new_status(i) = gen(i).QPform.Stor.UsableSize;
            elseif new_status(i)<0
                %Changed locked, but there was spare capacity in other timesteps so storage will not actually deplete
%                     disp(strcat('Warning, Potentially Depleating Storage in ',num2str(sum(dt(1:t))),'_hours'))
                new_status(i) = 0;
            end
        end
    end
    status = new_status;
    gen_output(t+1,:) = status;
end
end%Ends function update_storage

function [binary_comb,disp_comb,cost_comb,verified,i_best,locked,must_replace] = avoid_generation(qp,gen,date,i_best,k,t1,t2,start_cost,binary_comb,disp_comb,cost_comb,verified,gen_output,locked,dt)
%need to re-do this for variable time steps
%looking for enough slack in other generators to avoid the overdepleating storage at the minimum, and to get to the same final state
n_s = length(dt);
n_g = length(gen);
must_replace = true;
cap = gen(k).Output.Capacity*gen(k).Size;
[useful_stored_energy,~,stor] = stor_state(gen,gen_output,dt);
out = [];
if strcmp(gen(k).Type,'Chiller')
    out = 'C';
    eff = gen(k).Output.Cooling;
elseif strcmp(gen(k).Type,'Heater')
    out = 'H';
    eff = gen(k).Output.Heat;
elseif ismember(gen(k).Type,{'CHP Generator';'Electric Generator';})
    if isfield(gen(k).Output,'Electricity')
        eff = gen(k).Output.Electricity;
    else
        eff = gen(k).Output.DirectCurrent;
    end
    if isfield(stor,'DC')
        out = 'DC';
    elseif isfield(stor,'E')
        out = 'E';
    end
end
scale_cost = update_cost(date,gen);
[~,~,spare_gen_cumulative] = gen_limit(gen,gen_output,locked,dt);
rmv_cost = 0;
if isfield(stor,out) && any(useful_stored_energy.(out))>0
    margin_cost = marginal_cost(gen,gen_output,date);
    rem_gen = zeros(n_s,1);
    rem_heat = zeros(n_s,1);
    for t = t1:1:t2-1
        rem_gen(t:end) = rem_gen(t:end) + gen_output(t+1,k)*dt(t);%need to replace this energy in the storage by the end, and replace enough early so that UsefulStoredEnergy - remStor + makeup does not go negative
        rmv_cost = rmv_cost + scale_cost(t+1,k)*gen_output(t+1,k)./interp1(cap,eff,gen_output(t+1,k))*dt(t);
        if isfield(gen(k).QPform,'constCost')
            rmv_cost = rmv_cost + gen(k).QPform.constCost*scale_cost(t+1,k)*dt(t);
        end
        if strcmp(gen(k).Type,'Chiller') && isfield(gen(k).QPform.constDemand,'E')
            rmv_cost = rmv_cost + (gen(k).QPform.constDemand.E + gen_output(t+1,k)./interp1(cap,eff,gen_output(t+1,k)))*min(nonzeros(margin_cost.E.Cost.SpinReserve(:,t,1)))*dt(t);%%total electricity used by chiller * marginal cost
        end
        if strcmp(gen(k).Type,'CHP Generator')
            rem_heat(t:end) = rem_heat(t:end) + gen_output(t+1,k)*dt(t);
            cum_out = 0;
            j = 1;
            if gen_output(t+1,k)>0
                rem_heat(t:end) = rem_heat(t:end) - gen(k).QPform.constDemand.H;
            end
            states = gen(k).QPform.states(1:nnz(~cellfun('isempty',gen(k).QPform.states(:,end))),end);
            while j<=length(states)
                d = min(gen_output(t+1,k) - cum_out,gen(k).QPform.(states{j}).ub(2));
                cum_out = cum_out + d;
                rem_heat(t:end) = rem_heat(t:end) + d*gen(k).QPform.output.H(j,2);
                j = j+1;
            end
        end
    end
    useful = useful_stored_energy.(out);
    if strcmp(gen(k).Type,'CHP Generator') && ~isfield(stor,'H')
        useful2 = 0;
    elseif strcmp(gen(k).Type,'CHP Generator')
        useful2 = useful_stored_energy.H;
    end
    possible_to_avoid = false;
    new_stor_power = zeros(n_s,1);
    if (spare_gen_cumulative.(out)(end)+0.5*useful(end))>=rem_gen(end) && all(rem_gen-spare_gen_cumulative.(out)<useful)%it can reach at least half the planned final SOC, and it never lets the SOC dip below 0
        if ~strcmp(gen(k).Type,'CHP Generator') || (all(rem_heat-spare_gen_cumulative.H<useful2) && spare_gen_cumulative.H(end)>=rem_heat(end))
            possible_to_avoid = true;
            must_replace = false;
            for t = t1:1:t2-1
                locked(t+1,k) = false; %turn off
                new_stor_power(t) = disp_comb{t}(i_best(t),k);
                new_case = nonzeros((1:length(binary_comb{t}(:,1)))'.*ismember(binary_comb{t},locked(t+1,:),'rows'));
                if~isempty(new_case)
                    i_best(t) = new_case(1);
                    if ~verified{t}(new_case(1))
                        [cost_comb{t},verified{t},disp_comb{t}] = test_more_cases(qp{t},i_best(t),cost_comb{t},verified{t},disp_comb{t},binary_comb{t});
                    end
                else
                    disp_comb{t}(end+1,:) = disp_comb{t}(i_best(t),:);
                    binary_comb{t}(end+1,:) = binary_comb{t}(i_best(t),:);
                    cost_comb{t}(end+1) = cost_comb{t}(i_best(t))-1e-3; %make cheaper than standard
                    verified{t}(end+1) = false;
                    i_best(t) = length(verified{t});
                    disp_comb{t}(i_best(t),k) = 0;%set output of this generator to zero
                    binary_comb{t}(i_best(t),k) = 0;%set output of this generator to zero
                end
            end
        end
    end
    %%try and replace the generation with the cheapest generation from any other spare capacity (otherwise storage at end of horizon drops)
    if possible_to_avoid
        if any(rem_gen>useful)%Need to make up generation before storage would go negative
            t_limit = min(nonzeros((1:n_s)'.*(rem_gen>useful)));
        else
            t_limit = n_s;
        end
        make_up_gen = 0;
        rem_e = max(rem_gen);%removed energy (maximum = sum of generator output that is turned off * dt)
        sort_margin_cost = sort_mc(margin_cost,out,t_limit,k,dt);
        if ~isempty(sort_margin_cost) && sort_margin_cost(end,1)>=rem_e
            make_up_cost = interp1(sort_margin_cost(:,1),sort_margin_cost(:,2),rem_e);
            if (make_up_cost-rmv_cost)<start_cost(k)%if spareGen is less than StartCost
                r = 1;
                while make_up_gen<rem_e && r<length(sort_margin_cost(:,1))
                    t = sort_margin_cost(r,4);
                    add_gen_energy = min(sort_margin_cost(r,5)*dt(t),rem_e-make_up_gen);
                    disp_comb{t}(i_best(t),sort_margin_cost(r,3)) = disp_comb{t}(i_best(t),sort_margin_cost(r,3)) + add_gen_energy/dt(t);
                    %%edit storage dispatch at this timestep so that UpdateStorageState works correctly
                    disp_comb{t}(i_best(t),:) = stor_add(gen,disp_comb{t}(i_best(t),:),add_gen_energy,k);
                    make_up_gen = make_up_gen + add_gen_energy;
                    r = r+1;
                end
            end
        end
    end
end
end%Ends function avoid_generation

function [binary_comb,disp_comb,cost_comb,verified,i_best,locked,cant_move] = move_generation(qp,gen,i_best,k,t1,t2,binary_comb,disp_comb,cost_comb,verified,gen_output,locked,inc,dt)
%need to re-do this for variable time steps
[useful_stored_energy,~,~,spare_stor_cap,~] = stor_state(gen,gen_output,dt);
switch gen(k).Type
    case {'CHP Generator';'Electric Generator';}
        if isfield(useful_stored_energy,'DC')
            out = 'DC';
        elseif isfield(useful_stored_energy,'E')
            out = 'E';
        end
    case 'Chiller'
        out = 'C';
    case 'Heater'
        out = 'H';
    case 'Cooling Tower'
        out = 'CW';      
end
n_g = length(gen);
n_s = length(disp_comb);
cant_move = true;
i_best_alt = i_best;

if isfield(useful_stored_energy,out) && any(useful_stored_energy.(out))>0
    add_stor = zeros(n_s,1);
    rem_stor = zeros(n_s,1);
    stops = nonzeros((1:n_s)'.*(~locked(2:end,k) & locked(1:n_s,k)));
    n = t2-t1; %number of steps generator is on for
    if any(stops<t1)
        %it was on previously, try to move earlier if storage permits
        t_stop = stops(find(stops<t1,1,'last'));
        %check if there is a feasible option to leave this generator on at earlier time steps
        for j = 0:1:n-1
            opt = binary_comb{t_stop+j}; %all feasible options tested at this time
            locked_now = opt(i_best_alt(t_stop+j),:);
            locked_now(k) = true;
            new_locked = nonzeros((1:length(opt(:,1)))'.*ismember(opt,locked_now,'rows'));
            if~isempty(new_locked)
                i_best_alt(t_stop+j) = new_locked(1);
                if ~verified{t_stop+j}(new_locked(1))
                    [cost_comb{t_stop+j},verified{t_stop+j},disp_comb{t_stop+j}] = test_more_cases(qp{t_stop+j},i_best_alt(t_stop+j),cost_comb{t_stop+j},verified{t_stop+j},disp_comb{t_stop+j},opt);
                end
            end
            add_stor(t_stop+j:t1+j-1) = add_stor(t_stop+j:t1+j-1) + gen_output(t1+j+1,k);%need to hold this shifted energy in the storage from the previous time the generator shut down, until it had come on before
        end
        if ~any(i_best_alt(t_stop:t_stop+n-1) == i_best(t_stop:t_stop+n-1)) && all(add_stor<spare_stor_cap.(out))%possible to have generator on earlier
            for j = 0:1:n-1
                locked(t_stop+j+1,inc) = binary_comb{t_stop+j}(i_best_alt(t_stop+j),inc); %best alternative option tested at this time
            end
            for t = t_stop+n:t2
                locked(t+1,k) = false; %turn off
                new_case = nonzeros((1:length(binary_comb{t}(:,1)))'.*ismember(binary_comb{t},locked(t+1,:),'rows'));
                if~isempty(new_case)
                    i_best_alt(t) = new_case(1);
                    if ~verified{t}(new_case(1))
                        [cost_comb{t},verified{t},disp_comb{t}] = test_more_cases(qp{t},i_best_alt(t),cost_comb{t},verified{t},disp_comb{t},binary_comb{t});
                    end
                else
                    disp_comb{t}(end+1,:) = disp_comb{t}(i_best_alt(t),:);
                    binary_comb{t}(end+1,:) = binary_comb{t}(i_best_alt(t),:);
                    cost_comb{t}(end+1) = cost_comb{t}(i_best_alt(t))-1e-3; %make cheaper than standard
                    verified{t}(end+1) = false;
                    i_best_alt(t) = length(verified{t});
                    disp_comb{t}(i_best_alt(t),k) = 0;%set output of this generator to zero
                    binary_comb{t}(i_best_alt(t),k) = 0;%set output of this generator to zero
                end
            end
            i_best = i_best_alt; %use the alternative index
            cant_move = false;
        end
    end
    %if it cant move earlier, because of lack of storage space, try moving later
    if cant_move
        starts = nonzeros((1:n_s)'.*(locked(2:end,k) & ~locked(1:n_s,k)));
        t_start = starts(find(starts>t2,1,'first'))-n;%new starting point, n steps before its scheduled start
        if isempty(t_start)
            t_start = n_s-n+1;%push as late as possible
        end
        %check if there is a feasible option to leave this generator on at earlier time steps
        for j = 0:1:n-1
            opt = binary_comb{t_start+j}; %all feasible options tested at this time
            locked_now = opt(i_best_alt(t_start+j),:);
            locked_now(k) = true;
            new_locked = nonzeros((1:length(opt(:,1)))'.*ismember(opt,locked_now,'rows'));
            if~isempty(new_locked)
                i_best_alt(t_start+j) = new_locked(1);
            end
            rem_stor(t1+j:t_start+j-1) = rem_stor(t1+j:t_start+j-1) + gen_output(t1+j+1,k);%need to hold this shifted energy in the storage from the previous time the generator shut down, until it had come on before
        end
        if ~any(i_best_alt(t_start:t_start+n-1) == i_best(t_start:t_start+n-1)) && all(rem_stor>useful_stored_energy.(out))%possible to have generator on later
            i_best = i_best_alt; %use the alternative index
            for j = 0:1:n-1
                locked(t_start+j+1,inc) = binary_comb{t_start+j}(i_best(t_start+j),inc); %best alternative option tested at this time
                disp_comb{t_start+j}(i_best(t_start+j),k) = gen_output(t1+j+1,k);
            end
            for t = t1:t_start-1
                locked(t+1,k) = false; %turn 
                disp_comb{t}(i_best(t),k) = 0;%set output of this generator to zero
            end
            cant_move = false;
        end
    end
end
end%Ends function move_generation

function [i_best,locked,no_alt,disp_comb,cost_comb,verified] = alt_generation(qp,gen,i_best,k,t1,t2,start_cost,binary_comb,disp_comb,cost_comb,verified,gen_output,locked,dt)
[useful_stored_energy,stor_gen_avail,~] = stor_state(gen,gen_output,dt);
no_alt = true;
n_s = length(dt);

%determine any posible substitute generators
n_g = length(start_cost);
inc = false(1,n_g);
switch gen(k).Type
    case 'Chiller'
        out = 'C';
        for i = 1:1:n_g
            if strcmp(gen(i).Type,'Chiller') || (strcmp(gen(i).Type,'Thermal Storage') && strcmp(gen(i).Source,'Cooling'))
                inc(i) = true;
            end
        end
    case 'Cooling Tower'
        out = 'CW';
        for i = 1:1:n_g
            if strcmp(gen(i).Type,'Cooling Tower')
                inc(i) = true;
            end
        end
    case 'Heater'
        out = 'H';
        for i = 1:1:n_g
            if strcmp(gen(i).Type,'Heater') || strcmp(gen(i).Type,'CHP Generator')  || (strcmp(gen(i).Type,'Thermal Storage') && strcmp(gen(i).Source,'Heat'))
                inc(i) = true;
            end
        end
    case {'CHP Generator';'Electric Generator';}
        if isfield(useful_stored_energy,'DC')
            out = 'DC';
        elseif isfield(useful_stored_energy,'E')
            out = 'E';
        end
        for i = 1:1:n_g
            if strcmp(gen(i).Type,'CHP Generator') || strcmp(gen(i).Type,'Electric Generator') || strcmp(gen(i).Type,'Electric Storage')
                inc(i) = true;
            end
            if strcmp(gen(k).Type,'CHP Generator') && (strcmp(gen(i).Type,'Heater') || (strcmp(gen(i).Type,'Thermal Storage') && strcmp(gen(i).Source,'Heat')))
                inc(i) = true;
            end
        end
end
i_best_alt = i_best;
d_cost = zeros(t2-t1,1);
d_gen = zeros(t2-t1,1);
alt_locked = locked;
%for the time that the generator in question is on, try combinations that don't have it on.
for t = t1:1:t2-1
    opt = binary_comb{t}; %all options at this time
    %%run all unverified cases that replace the generator in question with either 1 or 2 or 3 others that are
    % a) running within (t2-t1) steps before/after and have a marginal cost less than starting the generator in question 
    % b) have a marginal cost + start-up cost less than starting the generator in question
    need_test = [];
    for c = 1:1:length(opt(:,1))
        if ~verified{t}(c) && ~opt(c,k)
            new_on = nonzeros((1:n_g).*((binary_comb{t}(c,:)-binary_comb{t}(i_best(t),:))>0));
            d_c = cost_comb{t}(c)-cost_comb{t}(i_best(t));
            if length(new_on)<=3 && all(inc(new_on)>0)
                for j = 1:1:length(new_on)
                    if ~any(locked(:,new_on(j))) 
                        d_c = d_c + start_cost(new_on(j));
                    end
                end
                if d_c<start_cost(k)
                    need_test(end+1) = c;
                end
            end
        end
    end
    [cost_comb{t},verified{t},disp_comb{t}] = test_more_cases(qp{t},need_test,cost_comb{t},verified{t},disp_comb{t},opt);
    
    %%determine the marginal cost of switching to the best possible option that doesn't use the generator in question
    cost2 = cost_comb{t} - cost_comb{t}(i_best(t)); %cost for these options
    cost2(opt(:,k)) = inf;%make the cost infinite for all combinations with generator i
    cost2(~verified{t}) = inf;%make the cost infinite for untested combinations
%     cost2(ismember(opt(:,inc),locked(t+1,inc),'rows')) = inf;%make the cost infinite if its not changing the status of any of the included generator type
    for j = 1:1:n_g
        if ~any(locked(t1:t2+1,j))%would be adding a startup
            cost2 = cost2 + opt(:,j)*start_cost(j)/(t2-t1);
        end
    end
    if any(~isinf(cost2))
        [d_cost(t-t1+1),i_best_alt(t)] = min(cost2);
        d_gen(t-t1+1) = disp_comb{t}(i_best(t),k) - disp_comb{t}(i_best_alt(t),k);
        alt_locked(t+1,:) = binary_comb{t}(i_best_alt(t),:);
    else
        d_cost = inf;%no feasible combinations without turning on a different generator (don't make changes to this segment)
        break
    end
end

%%If it wants to swap to have a different generator on now, don't constrain its maximum output by shutting it down immediately
new_on = zeros(n_s,n_g);
for t = t1:1:t2-1
    new_on(t,:)= ((binary_comb{t}(i_best_alt(t),:)-binary_comb{t}(i_best(t),:))>0);
end
for i = 1:1:n_g
    if all(new_on(t1:t2-1,i))
        alt_locked(t2:end,i) = true;
    end
end
%%check if the alternate operation strategy is feasible when put in sequence
spare_gen = zeros(t2-t1,1);
max_out = gen_limit(gen,gen_output,alt_locked,dt);
heat_out = zeros(t2-t1,n_g);
for i = 1:1:n_g
    if strcmp(gen(i).Type,'CHP Generator')
        heat_out(:,i) = chp_heat(gen(i).QPform,gen_output(t1+1:t2,i));%convert GenOutput to heatOut
    end
    if i~=k && disp_comb{t}(i_best_alt(t),i)>0
        if strcmp(gen(k).Type,gen(i).Type) || (strcmp(gen(k).Type,'CHP Generator') && strcmp(gen(i).Type,'Electric Generator')) || (strcmp(gen(k).Type,'Electric Generator') && strcmp(gen(i).Type,'CHP Generator'))
            spare_gen = spare_gen + max_out.(out)(t1+1:t2,i) - gen_output(t1+1:t2,i);%capacity is either UB or limited by ramping
        end
    end
    if strcmp(gen(k).Type,'Heater') && strcmp(gen(i).Type,'CHP Generator')
        for t = t1:1:t2-1
            spare_gen(t-t1+1) = spare_gen(t-t1+1) + max_out.H(t+1,i) - heat_out(t-t1+1,i);
        end
    end  
end
if strcmp(gen(k).Type,'CHP Generator')
    d_heat = heat_out(:,k);
    spare_heat = zeros(t2-t1,1);
    heat_out(:,k) = 0;
    for i = 1:1:n_g
        if strcmp(gen(i).Type,'Heater')
            spare_heat= spare_heat + max_out.H(t1+1:t2,i) - gen_output(t1+1:t2,i);
        elseif strcmp(gen(i).Type,'CHP Generator')
            spare_heat= spare_heat + max_out.H(t1+1:t2,i) - heat_out(:,i);
        end
    end
end
useful1 = 0;
useful2 = 0;
if isfield(useful_stored_energy,out)
    useful1 = stor_gen_avail.(out)(t1:t2-1);
    useful2 = min(useful_stored_energy.(out)(t1:end));    
end
if sum(d_cost)<start_cost(k) && all(d_gen<(spare_gen+useful1./dt(t1:t2-1))) && sum(d_gen-spare_gen)<useful2 %sum of the marginal increase in cost is less than the start-up cost, there is spare capacity in the other generators & storage, and the cumulative loss of generation does not deplete the storage
    if strcmp(gen(k).Type,'CHP Generator')
        useful3 = 0;
        useful4 = 0;
        if isfield(useful_stored_energy,'H')
            useful3 = stor_gen_avail.H(t1:t2-1);
            useful4 = min(useful_stored_energy.H(t1:end));  
        end
    end
    if ~strcmp(gen(k).Type,'CHP Generator') || (all(d_heat<(spare_heat+useful3./dt(t1:t2-1))) && sum(d_heat-spare_heat)<useful4)
        i_best = i_best_alt; %use the alternative index
        for t = t1:1:t2-1
            locked(t+1,inc) = binary_comb{t}(i_best(t),inc); %best alternative option tested at this time
        end
        no_alt = false;
    end
end
end%Ends function alt_generation

function [i_best,locked,cant_keep_on,disp_comb,cost_comb,verified] = leave_gen_on(qp,gen,i_best,k,t1,t2,start_cost,binary_comb,disp_comb,cost_comb,verified,gen_output,locked,dt)
%% Find the cheapest feasible alternative dispatch that keeps this generator on (only use generators that were on previously or will be on)
[~,~,~,spare_stor_cap,stor_slack_avail] = stor_state(gen,gen_output,dt);
switch gen(k).Type
    case {'CHP Generator';'Electric Generator';}
        if isfield(spare_stor_cap,'DC')
            out = 'DC';
        elseif isfield(spare_stor_cap,'E')
            out = 'E';
        end
    case 'Chiller'
        out = 'C';
    case 'Heater'
        out = 'H';
    case 'Cooling Tower'
        out = 'CW';      
end
cant_keep_on = true;
%only allow generators that are on at begining or end to be involved (or that have smaller start-up cost)
i_best_alt = i_best;
d_cost = zeros(t2-t1,1);
d_gen = zeros(t2-t1,1);
slack_gen = zeros(t2-t1,1);
if t1 == 1
    prev = gen_output(1,:);
else
    prev = disp_comb{t1-1}(i_best(t1-1),:);
end
for t = t1:1:t2-1
    opt = binary_comb{t}; %all feasible options tested at this time
    %%run unverified cases that just add the generator in question
    locked_now = opt(i_best_alt(t),:);
    locked_now(k) = true;
    new_locked = nonzeros((1:length(opt(:,1)))'.*ismember(opt,locked_now,'rows'));
    if~isempty(new_locked)
        i_best_alt(t) = new_locked(1);
        if ~verified{t}(new_locked(1))
            [cost_comb{t},verified{t},disp_comb{t}] = test_more_cases(qp{t},i_best_alt(t),cost_comb{t},verified{t},disp_comb{t},opt);
        end
        d_cost(t-t1+1) = cost_comb{t}(i_best_alt(t)) - cost_comb{t}(i_best(t));
        [d_gen(t-t1+1),slack_gen(t-t1+1)] = slack_cap(gen,disp_comb{t}(i_best(t),:),disp_comb{t}(i_best_alt(t),:),prev,k,sum(dt(t1:t)));
    else
        d_cost = inf;%no feasible combinations keeping this generator active (don't make changes to this segment)
        break
    end
end
useful1 = 0;
useful2 = 0;
if isfield(spare_stor_cap,out)
    useful1 = stor_slack_avail.(out)(t1:t2-1);
    useful2 = min(spare_stor_cap.(out)(t1:end));
end
if sum(d_cost)<start_cost(k) && all(d_gen<(slack_gen+useful1./dt(t1:t2-1))) && sum(d_gen-slack_gen)<useful2 %sum of the marginal increase in cost is less than the start-up cost, there is spare capacity in the other generators & storage, and the cumulative loss of generation does not deplete the storage
    i_best = i_best_alt; %use the alternative index
    for t = t1:1:t2-1
        locked(t+1,:) = binary_comb{t}(i_best(t),:); %best alternative option tested at this time
    end
    cant_keep_on = false;
end
end%Ends function leave_gen_on

function [on_seg, off_seg] = seg_length(locked,start_cost,dt,skip_on,skip_off,inc)
on_seg = [];
off_seg = [];
n_g = length(start_cost);
n_s = length(dt);
% find length of segments that a generator is on or off
for i = 1:1:n_g
    if start_cost(i)>0 && any(~locked(:,i)) && inc(i)
        starts = nonzeros((1:n_s)'.*(locked(2:end,i) & ~locked(1:n_s,i)));
        stops = nonzeros((1:n_s)'.*(~locked(2:end,i) & locked(1:n_s,i)));
        n_on = length(starts);
        n_off = length(stops);
        if n_on>0 && n_off>0 %only look at generators that turn both on and off during the window
            if stops(1)<starts(1) %generator is off for a segment
                seg = [i, stops(1), starts(1), sum(dt(stops(1):starts(1)-1))];
                if (isempty(skip_off) || ~any(ismember(skip_off,seg,'rows')))% && Seg(4)<Horizon/4 %segment hasn't previously been ruled out, and is less than 1/4 of the horizon
                    off_seg(end+1,1:4) = seg;
                end
                stops = stops(2:end);
                n_off  = n_off - 1;
            end
            j = 0;
            while j < n_off
                j = j + 1;
                seg = [i, starts(j), stops(j), sum(dt(starts(j):stops(j)-1))]; %index of generator, start index, stop index, duration of segment
                if (isempty(skip_on) || ~any(ismember(skip_on,seg,'rows')))% && Seg(4)<Horizon/4 %segment hasn't previously been ruled out, and is less than 1/4 of the horizon
                    on_seg(end+1,1:4) = seg;
                end
                if j<n_on
                    seg = [i, stops(j), starts(j+1), sum(dt(stops(j):starts(j+1)-1))];
                    if (isempty(skip_off) || ~any(ismember(skip_off,seg,'rows')))% && Seg(4)<Horizon/4 %segment hasn't previously been ruled out, and is less than 1/4 of the horizon
                        off_seg(end+1,1:4) = seg;
                    end
                end
            end
        end
    end
end
end%Ends function seg_length

function [d_gen,slack_gen] = slack_cap(gen,original,new,prev,k,dt)
n_g = length(gen);
d_gen = gen(k).QPform.A.lb(end);
if strcmp(gen(k).Type,'Chiller')
    include = {'Chiller'};
elseif strcmp(gen(k).Type,'Heater')
    include = {'Heater'};
elseif ismember(gen(k).Type,{'CHP Generator';'Electric Generator';})
    include = {'CHP Generator';'Electric Generator';};
end
slack_gen = 0;
for i = 1:1:n_g
    if i~=k && ismember(gen(i).Type,include) 
        if original(i)>0 && new(i) >0
            slack_gen = slack_gen + min(original(i)-gen(i).QPform.A.lb(end),original(i)-prev(i)+gen(i).VariableStruct.dX_dt*dt);%slack capacity is either Original - LB or limited by ramping
        elseif original(i)>0
            slack_gen = slack_gen + original(i);
        elseif new(i)>0
            slack_gen = slack_gen - new(i);
        end
    end
end
end%Ends function slack_cap

function [stored_energy,stor_gen_avail,stor,spare_stor_cap,stor_slack_avail] = stor_state(gen,gen_output,dt)
n_g = length(gen);
n_s = length(gen_output(2:end,1));
stored_energy =[];
for i = 1:1:n_g
    if isfield(gen(i).QPform,'Stor') && ~strcmp(gen(i).Type,'Hydro Storage')
        out = char(fieldnames(gen(i).QPform.output));
        if ~isfield(stored_energy,out)
            stored_energy.(out) = zeros(n_s,1);
            stor_gen_avail.(out) = zeros(n_s,1);
            spare_stor_cap.(out) = zeros(n_s,1);
            stor_slack_avail.(out) = zeros(n_s,1);
            stor.(out) = [];
        end
        buff = gen(i).QPform.Stor.UsableSize*(gen(i).VariableStruct.Buffer/100);
        stored_energy.(out) = stored_energy.(out) + (gen_output(2:end,i)-buff)*gen(i).QPform.Stor.DischEff;
        stor_gen_avail.(out) = stor_gen_avail.(out) +  min((gen_output(2:end,i)-buff)*gen(i).QPform.Stor.DischEff./dt,gen(i).QPform.Stor.PeakDisch);
        spare_stor_cap.(out) = spare_stor_cap.(out) + (gen(i).QPform.Stor.UsableSize - buff - gen_output(2:end,i))/gen(i).QPform.Stor.ChargeEff;
        stor_slack_avail.(out) = stor_slack_avail.(out) +  min((gen(i).QPform.Stor.UsableSize-buff-gen_output(2:end,i))/gen(i).QPform.Stor.ChargeEff./dt,gen(i).QPform.Stor.PeakCharge);
        stor.(out)(end+1) = i;
    end
end
end%Ends function stor_state


function dispatch = stor_add(gen,prev,add_gen,k)
n_g = length(gen);
dispatch = prev;
for i = 1:1:n_g
    if strcmp(gen(k).Type,'Chiller') && strcmp(gen(i).Type,'Thermal Storage') && strcmp(gen(i).Source,'Cooling')
        change = min(add_gen, gen(i).QPform.Stor.PeakCharge + prev(i));
        add_gen = add_gen - change;
        dispatch(i) = prev(i) - change;
    end
    if strcmp(gen(k).Type,'Heater') && strcmp(gen(i).Type,'Thermal Storage') && strcmp(gen(i).Source,'Heat')
        change = min(add_gen, gen(i).QPform.Stor.PeakCharge + prev(i));
        add_gen = add_gen - change;
        dispatch(i) = prev(i) - change;
    end
    if ismember(gen(k).Type,{'CHP Generator';'Electric Generator';}) && strcmp(gen(i).Type,'Electric Storage')
        change = min(add_gen, gen(i).QPform.Stor.PeakCharge + prev(i));
        add_gen = add_gen - change;
        dispatch(i) = prev(i) - change;
    end
end
end%Ends function stor_add

function [cost_comb,verified,disp_comb] = test_more_cases(qp,need_test,cost_comb,verified,disp_comb,opt)
parallel = false;
if license('test','Distrib_Computing_Toolbox') 
    parallel = true;
end
if ~isempty(need_test)
    n_test = length(need_test);
    if parallel %&& n_test>10
        cost_par = cost_comb(need_test);
        disp_comb_par = disp_comb(need_test,:);
        opt_par = opt(need_test,:);
        qp_par = qp;
        parfor par_i = 1:n_test
            qp_test_par = disable_generators_step(qp_par,opt_par(par_i,:)>0);%Disable generators here
            [x_par, flag1] = call_solver(qp_test_par);
            if flag1 == 1
                cost_par(par_i) = 0.5*x_par'*qp_test_par.H*x_par + x_par'*qp_test_par.f + sum(qp_test_par.constCost.*(opt_par(par_i,:)>0));
                disp_comb_par(par_i,:) = sort_solution_step(x_par,qp_test_par);
            else
                cost_par(par_i) = inf;
            end
        end
        cost_comb(need_test) = cost_par;
        verified(need_test) = true;
        disp_comb(need_test,:) = disp_comb_par;
    else
        for j = 1:1:length(need_test)
            c = need_test(j);
            qp_test = disable_generators_step(qp,opt(c,:)>0);%Disable generators here
            [x, flag1] = call_solver(qp_test);
            if flag1 == 1
                cost_comb(c) = 0.5*x'*qp_test.H*x + x'*qp_test.f + sum(qp_test.constCost.*(opt(c,:)>0));
                disp_comb(c,:) = sort_solution_step(x,qp_test);
            else
                cost_comb(c) = inf;
            end
        end
    end
end
end%ends function test_more_cases