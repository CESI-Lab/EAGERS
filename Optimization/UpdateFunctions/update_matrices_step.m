function qp = update_matrices_step(gen,building,fluid_loop,subnet,options,qp,forecast,scale_cost,marginal,stor_power,dt,ic,first_profile,limit,t,temperatures)
% update the equalities with the correct demand, and scale fuel and electric costs
% ec is the expected end condition at this time stamp (can be empty)
% stor_power is the expected output/input of any energy storage at this timestep (can be empty)
% min_power and MaxPower define the range of this generator at this timestep
% temperatures is the building and cooling tower water loop temperatures
qp.solver = options.solver;
n_g = length(gen);

%% determine which components are enabled
qp.Organize.Enabled = true(1,n_g);
for i = 1:1:n_g
    if ~gen(i).Enabled
        qp.Organize.Enabled(i) = false;
    end
end
%% update the demands and self-discharge losses
qp = update_mat_demand_step(gen,qp,subnet,options,forecast,first_profile,ic,dt,t);

%% Building inequalities & equality
qp = update_mat_building_step(building,qp,forecast,temperatures,dt,t);
%% Fluid loop equality
qp = update_fluid_loop_step(fluid_loop,qp,temperatures,dt,t);

%% Update upper and lower bounds based on ramping constraint (if applicable)
[min_power,max_power] = constrain_min_max(gen,limit,ic,dt,t);
qp = update_gen_limit(gen,qp,min_power,max_power);

%% update storage bounds
if ~isempty(first_profile)
    qp = update_storage_bounds(gen,qp,stor_power,first_profile,dt,t);
end

%% Adding cost for Line Transfer Penalties to differentiate from spillway flow 
qp = update_mat_hydro_step(gen,qp,subnet);

%% update costs
qp = update_mat_cost_step(gen,qp,first_profile,stor_power,scale_cost,marginal,dt,t);

%% update spinning reserve
qp = update_spin_reserve_step(qp,options,forecast,dt,t);

end%ends function update_matrices_step

function qp = update_mat_demand_step(gen,qp,subnet,options,forecast,first_profile,ic,dt,t)
%% update demands and storage self-discharge
n_g = length(gen);
network_names = fieldnames(subnet);
for nn = 1:1:length(network_names)
    net = network_names{nn};
    abrev = subnet.(net).abbreviation;
    qp.network.(abrev).nodes = length(subnet.(net).nodes);
    if strcmp(net,'Hydro')
        %don't do a water balance, since it depends on multiple time steps.
        %Any extra outflow at this time step is subtracted from the expected
        %outflow at the next step (same SOC and flows up river).
    else
        if strcmp(net,'Electrical') && options.SpinReserve
            qp.b(qp.Organize.SpinReserve) = -options.SpinReservePerc/100*sum(forecast.Demand.E(t,:),2);% -shortfall + SRancillary - SR generators - SR storage <= -SR target
        end
        qp.constDemand.(abrev).req = zeros(qp.network.(abrev).nodes,1);
        qp.constDemand.(abrev).load = zeros(qp.network.(abrev).nodes,n_g);
        for n = 1:1:qp.network.(abrev).nodes
            equip = subnet.(net).Equipment{n};
            req = qp.Organize.Balance.(net)(n);
            loads = nonzeros(qp.Organize.Demand.(net){n});
            if ~isempty(loads)%need this to avoid the next line when there is no load
                qp.beq(req) = sum(forecast.Demand.(abrev)(t,loads),2); %multiple demands can be at the same node
            end
            for j = 1:1:length(equip)
                k = equip(j);
                if strcmp(net,'Electrical') 
                    if strcmp(gen(k).Source,'Renewable')% subtract renewable generation 
                        qp.beq(req) = qp.beq(req) - forecast.Renewable(t,k); %put renewable generation into energy balance at correct node
                    end
                end
                if ~isempty(first_profile) && ismember(gen(k).Type,{'Electric Storage';'Thermal Storage';})
                    if isfield(gen(k).QPform.output,abrev)
                        loss = (gen(k).QPform.Stor.SelfDischarge*gen(k).QPform.Stor.UsableSize*gen(k).QPform.Stor.DischEff);
                        d_soc = first_profile(t+1,k) - ic(k);%positive means charging
                        qp.beq(req) = qp.beq(req) + (d_soc/dt(t)+loss)*gen(k).QPform.Stor.DischEff;
                        qp.b(qp.Organize.Inequalities(k)) = -(1/gen(k).QPform.Stor.ChargeEff-gen(k).QPform.Stor.DischEff)*(d_soc/dt(t)+loss);
                    end
                end
                if isfield(gen(k).QPform,'constDemand') && isfield(gen(k).QPform.constDemand,abrev)
                    qp.constDemand.(abrev).req(n,1) = req;
                    qp.constDemand.(abrev).load(n,k) = gen(k).QPform.constDemand.(abrev);
                end
            end
        end
    end
end  
end%Ends function update_mat_demand_step

function qp = update_mat_cost_step(gen,qp,first_profile,stor_power,scale_cost,marginal,dt,t)
n_g = length(gen);
H = diag(qp.H);%convert to vector
for i = 1:1:n_g
    states = qp.Organize.States{i}; %states of this generator
    switch gen(i).Type
        case 'Utility'
            if strcmp(gen(i).Source,'Electricity')
                qp.f(states(1)) = qp.f(states(1))*scale_cost(i)*dt(t);
                if gen(i).VariableStruct.SellBackRate == -1
                    qp.f(states(2)) = qp.f(states(2))*scale_cost(i)*dt(t); %sellback is a fixed percent of purchase costs (percent set wehn building matrices)
                else
                    %constant sellback rate was taken care of when building the matrices
                end
            end
        case {'Electric Generator';'CHP Generator';'Chiller';'Heater';}%all generators and utilities
            H(states) = H(states)*scale_cost(i)*dt(t);
            qp.f(states) = qp.f(states)*scale_cost(i)*dt(t);
        case {'Electric Storage';'Thermal Storage';'Hydro Storage';}%all  storage
            stor_cat = fieldnames(gen(i).QPform.output);
            %penalize deviations from the expected storage output
            a = 1; %sets severity of quadratic penalty
            %scale the penalty by a) larger of 5% capcity and 2x the
            %expected change in stored capacity, b) lesser of a) and
            %remaining space in storage
            MaxCharge = 0.05*gen(i).QPform.Stor.UsableSize/dt(t);
            if ~isempty(first_profile)
                MaxCharge = max(MaxCharge,2*abs(stor_power(i)));
                if gen(i).QPform.Stor.UsableSize>first_profile(t+1,i)
                    MaxCharge = min((gen(i).QPform.Stor.UsableSize-first_profile(t+1,i))/dt(t),MaxCharge);
                end
                if first_profile(t+1,i)>0
                    MaxCharge = min(0.5*first_profile(t+1,i)/dt(t),MaxCharge);
                end
            end
            qp.f(states(1)) = marginal.(stor_cat{1});
            H(states(1)) = a*2*qp.f(states(1))/MaxCharge;  %factor of 2 because its solving C = 0.5*x'*H*x + f'*x
            if length(states)>1
                qp.f(states(2)) = 1e-6; %small penalty on charging state so that if there are no other generators it keeps the charging constraint as an equality
            end
    end 
end 
qp.H = diag(H);%convert back to matrix
qp.constCost = qp.constCost.*scale_cost*dt(t);
end%Ends function update_mat_cost_step

function qp = update_spin_reserve_step(qp,options,forecast,dt,t)
n_g = length(qp.Organize.Dispatchable);
if options.SpinReserve
    sr_short = qp.Organize.SpinReserveStates(1,n_g+1);%cumulative spinning reserve shortfall 
    sr_ancillary = qp.Organize.SpinReserveStates(1,n_g+2);%Ancillary spinning reserve value (negative cost) 
    if options.SpinReservePerc>5 %more than 5% spinning reserve
        spin_cost = 2*dt(t)./(options.SpinReservePerc/100*sum(forecast.Demand.E(t,:),2));% -shortfall + SRancillary - SR generators - SR storage <= -SR target
    else
        spin_cost = 2*dt(t)./(0.05*sum(forecast.Demand.E(t,:),2));% -shortfall + SRancillary - SR generators - SR storage <= -SR target
    end
    qp.H(sr_short,sr_short) = spin_cost;%effectively $2 per kWh at point where shortfall = spin reserve percent*demand or $2 at 5%
    qp.f(sr_short) = 0.05*dt(t); % $0.05 per kWh
end
end%Ends function update_spin_reserve_step

function qp = update_mat_building_step(building,qp,forecast,temperatures,dt,t)
n_b = length(building);
n_g = length(qp.Organize.Dispatchable);
n_l = length(qp.Organize.Transmission);
for i = 1:1:n_b
    states = qp.Organize.States{n_g+n_l+i};
    %Heating cooling offsets (value through dead_band)    
    qp.Organize.Building.H_Offset(1,i) = forecast.Building.H0(t,i) - max(0,building(i).QPform.UA*(forecast.Building.Tzone(t,i)-forecast.Building.Tset_H(t,i)));
    qp.Organize.Building.C_Offset(1,i) = forecast.Building.C0(t,i) - max(0,building(i).QPform.UA*(forecast.Building.Tset_C(t,i)-forecast.Building.Tzone(t,i)));
    
    %electrical energy balance
    req = qp.Organize.Building.Electrical.req(i);
    qp.beq(req) = qp.beq(req) + forecast.Building.E0(t,i);%Equipment and nominal Fan Power
    %ensure when heating = H0, that there is no change in electrical demand
    qp.beq(req) = qp.beq(req) + (forecast.Building.H0(t,i) - qp.Organize.Building.H_Offset(1,i))*qp.Aeq(req,states(2));
    %ensure when cooling = C0, that there is no change in electrical demand
    qp.beq(req) = qp.beq(req) + (forecast.Building.C0(t,i) - qp.Organize.Building.C_Offset(1,i))*qp.Aeq(req,states(3));
    
    %Heating Equality
    req = qp.Organize.Building.DistrictHeat.req(i);
    qp.beq(req) = qp.beq(req) + qp.Organize.Building.H_Offset(1,i);%Heating load
%     qp.lb(states(2)) = qp.lb(states(2)) + qp.Organize.Building.H_Offset(1,i);
    
    %Cooling Equality
    req = qp.Organize.Building.DistrictCool.req(i);
    qp.beq(req) = qp.beq(req) + qp.Organize.Building.C_Offset(1,i);%Cooling Load
%     qp.lb(states(3)) = qp.lb(states(3)) + qp.Organize.Building.C_Offset(1,i);
    
    %heating inequality @ T0  H_bar >= H0 - offset becomes: 
    qp.b(qp.Organize.Building.r(i)) = (building(i).QPform.UA + building(i).QPform.Cap/(3600*dt(t)))*temperatures.build(1,i) - (forecast.Building.H0(t,i) - qp.Organize.Building.H_Offset(1,i));
    qp.A(qp.Organize.Building.r(i),states(1)) = (building(i).QPform.UA+building(i).QPform.Cap/(3600*dt(t)));
    %cooling inequality @ T0  C_bar >= C0 - offset
    qp.b(qp.Organize.Building.r(i)+1) = (building(i).QPform.UA + building(i).QPform.Cap/(3600*dt(t)))*temperatures.build(1,i) - (forecast.Building.C0(t,i) - qp.Organize.Building.C_Offset(1,i));
    qp.A(qp.Organize.Building.r(i)+1,states(1)) = -building(i).QPform.UA - building(i).QPform.Cap/(3600*dt(t));
    %penalty states (excess and under temperature)
    qp.b(qp.Organize.Building.r(i)+2) = forecast.Building.Tmax(t,i);%upper buffer inequality, the longer the time step the larger the penaly cost on the state is.
    qp.b(qp.Organize.Building.r(i)+3) = -forecast.Building.Tmin(t,i);%lower buffer inequality
end
end%Ends function update_mat_building_step

function qp = update_fluid_loop_step(fluid_loop,qp,temperatures,dt,t)
n_b = length(qp.Organize.Building.r);
n_g = length(qp.Organize.Dispatchable);
n_l = length(qp.Organize.Transmission);
n_fl = length(fluid_loop);
for i = 1:1:n_fl
    %1st order temperature model: 0 = -T(k) + T(k-1) + dt/Capacitance*(energy balance)
    state = qp.Organize.States{n_g+n_l+n_b+i};
    req = qp.Organize.Balance.CoolingWater(i);
    qp.Organize.fluid_loop(i) = req;
    capacitance = fluid_loop(i).fluid_capacity*fluid_loop(i).fluid_capacitance; %Water capacity in kg and thermal capacitance in kJ/kg*K to get kJ/K
    qp.Aeq(req,:) = qp.Aeq(req,:); %energy balance for chillers & cooling tower fans already put in this row of Aeq
    qp.Aeq(req,state) = -capacitance/(3600*dt(t)); %-T(k)
    qp.beq(req) = -temperatures.fluid_loop(i)*capacitance/(3600*dt(t));%-T(k-1)
end
end%Ends function update_fluid_loop_step


function [min_power,max_power] = constrain_min_max(gen,limit,ic,dt,t)
n_g = length(gen);
upper_bound = zeros(1,n_g);
dx_dt = zeros(1,n_g);
for i = 1:1:n_g
    switch gen(i).Type
        case 'Hydro Storage'
            %%do something to set UB
        case 'Utility'
            if strcmp(gen(i).Source,'Electricity')
                states = gen(i).QPform.states;
                for j = 1:1:length(states)
                    upper_bound(i) = upper_bound(i) + gen(i).QPform.(states{j}).ub(end);
                end
            end
        case {'Electric Generator';'CHP Generator';'Chiller';'Heater';}
            states = gen(i).QPform.states(1:nnz(~cellfun('isempty',gen(i).QPform.states(:,end))),end);
            for j = 1:1:length(states)
                upper_bound(i) = upper_bound(i) + gen(i).QPform.(states{j}).ub(end);
            end
        otherwise
            %do nothing
    end
    if isfield(gen(i).VariableStruct,'dX_dt')
        dx_dt(i) = gen(i).VariableStruct.dX_dt;
    end
end
if strcmp(limit, 'unconstrained')
    min_power = 0*upper_bound;
    max_power = upper_bound;
elseif strcmp(limit, 'constrained')
    min_power = max(0,ic(1:n_g)-dx_dt*dt(t));
    max_power = min(upper_bound,ic(1:n_g)+dx_dt*dt(t));
else
    min_power = max(0,ic(1:n_g)-dx_dt*sum(dt(1:t)));
    max_power = min(upper_bound,ic(1:n_g)+dx_dt*sum(dt(1:t)));
end
end%Ends function constrain_min_max

function qp = update_gen_limit(gen,qp,min_power,max_power)
n_g = length(gen);
for i = 1:1:n_g
    switch gen(i).Type
        case {'Utility';'Electric Generator';'CHP Generator';'Chiller';'Heater';}%all generators and utilities
            states = qp.Organize.States{i};
            if min_power(i)>sum(qp.lb(states))%raise the lower bounds
                qp.Organize.Dispatchable(i) = false; %can't shut off
                power_add = min_power(i) - sum(qp.lb(states));
                for f = 1:1:length(states) %starting with lb of 1st state in generator, then moving on
                    gap = qp.ub(states(f)) - qp.lb(states(f));
                    if gap>0
                        gap = min(gap,power_add);
                        qp.lb(states(f)) = qp.lb(states(f)) + gap;
                        power_add = power_add-gap;
                    end
                end
            end
            if max_power(i)<gen(i).Size %lower the upper bounds
                power_sub = sum(qp.ub(states)) - max_power(i);
                for f = length(states):-1:1 %starting with lb of 1st state in generator, then moving on
                    if qp.ub(states(f))>0
                        gap = min(qp.ub(states(f)),power_sub);
                        qp.ub(states(f)) = qp.ub(states(f)) - gap;
                        power_sub = power_sub-gap;
                    end
                end
            end
        case {'Hydro Storage';}
            states = qp.Organize.States{i};
            qp.lb(states(1)) = min_power(i);
            qp.ub(states(1)) = max_power(i);
    end
end
end%Ends function update_gen_limit

function qp = update_storage_bounds(gen,qp,stor_power,first_profile,dt,t)
n_g = length(gen);
for i = 1:1:n_g
    states = qp.Organize.States{i};
    switch gen(i).Type
        case {'Electric Storage';'Thermal Storage';'Hydrogen Storage';}
            %update storage output range to account for what is already scheduled
            f = fieldnames(gen(i).QPform.output);
            switch f{1}
                case 'H'
                    req = qp.Organize.Balance.DistrictHeat; %rows of Aeq associated with heat demand
                case 'E'
                    req = qp.Organize.Balance.Electrical;
                case 'DC'
                    req = qp.Organize.Balance.DirectCurrent;
                case 'C'
                    req = qp.Organize.Balance.DistrictCool; %rows of Aeq associated with cool demand
            end
            s = states(1);
            chargingSpace = max(0,sum((gen(i).QPform.Stor.UsableSize-first_profile(t+1,i)).*qp.Aeq(req,s))/dt(t));%you can't charge more than you have space for            
            qp.lb(s) = max(qp.lb(s) - stor_power(i), -chargingSpace);%change in storage for this power output
            qp.ub(s) = min(qp.ub(s) - stor_power(i), first_profile(t+1,i)/dt(t));%you can't discharge more than you have stored            
            if qp.Organize.SpinRow(i)~=0 %update spinning reserve (max additional power from storage
                qp.beq(qp.Organize.SpinRow{i}) = qp.ub(s);
            end
        case 'Hydro Storage'
            s = states(1);
            qp.lb(s) = qp.lb(s) - stor_power(i);%change in storage for this power output
            qp.ub(s) = qp.ub(s) - stor_power(i);%change in storage for this power output
            %update the equality with NewPower*conversion + spill - Outflow = nominal PowerGen Flow
            h = nonzeros((1:length(qp.Organize.Balance.Hydro))'.*(i==qp.Organize.Balance.Hydro));
            qp.beq(qp.Organize.HydroEqualities(h)) = -first_profile(i)*gen(i).QPform.Stor.Power2Flow;
    end
    qp.lb(states) = min(qp.lb(states),qp.ub(states)); %fix rounding errors
end
end%Ends function update_storage_bounds

function qp = update_mat_hydro_step(gen,qp,subnet)
%%also want to penalize spill flow to conserve water for later, spill flow
%%is just to meet instream requirments when the power gen needs to be low
if isfield(subnet,'Hydro') && ~isempty(subnet.Electrical.lineNames)
    n_g = length(gen);
    for k = 1:1:length(subnet.Electrical.lineNames)
        line = subnet.Electrical.lineNumber(k);
        states = qp.Organize.States{n_g+line};
        if length(states)>1 %bi-directional transfer with penalty states
            s2 = states(2):qp.Organize.t1States:(states(2));%penalty from a to b
            s3 = states(3):qp.Organize.t1States:(states(3));% penalty from b to a
            qp.f(s2) = 0.001; %Adding .1 cent/kWhr cost to every Electric line penalty to differentiate from spill flow
            qp.f(s3) = 0.001; %Adding .1 cent/kWhr cost to every Electric line penalty to differentiate from spill flow
        end
    end 
    for i = 1:1:n_g
        if strcmp(gen(i).Type,{'Hydro Storage';}) && isfield(gen(i).QPform,'S') 
            states = qp.Organize.States{i};
            qp.f(states(2)) = 0.0001/gen(i).QPform.Stor.Power2Flow;
        end
    end
end 
end%Ends function update_mat_hydro_step