function qp = update_matrices(gen,buildings,fluid_loop,subnet,options,qp,date,scale_cost,margin_cost,forecast,ec)
% gen is the set of generators in the project
% building is the set of buildings in the project
% fluid_loop is the set of fluid loops, e.g. cooling tower water loop
% subnet is the agregated networks
% options are the optimization options
% qp is the optimization problem 
% date is the vector of time steps
% scale_cost is the matrix multiplier of cost for each generator at each time step
% margin_cost is the marginal cost of generation (used in the final value of the energy storage)
% forecast is the forecasted loads
% ec is the end condition for the threshold optimization
if strcmp(options.solver,'ANN')
    qp.solver = 'quadprog';
else
    qp.solver = options.solver;
end
n_g = length(gen);
n_b = length(buildings);

dt = (date(2:end) - date(1:end-1))*24;
n_s = length(dt);
n_l = length(qp.Organize.Transmission);
n_fl = length(fluid_loop);

%update initial conditions
qp.Organize.Enabled = true(1,n_g);
for i = 1:1:n_g
    if qp.Organize.IC(i)>0 %generators and storage devices
        qp.beq(qp.Organize.IC(i)) = gen(i).CurrentState(1);
        qp.ub(qp.Organize.IC(i)) = gen(i).CurrentState(1)+1;%+1 just to help solver find feasible
    end
    if ~gen(i).Enabled
        qp.Organize.Enabled(i) = false;
    end
end
for i = 1:1:n_b
    if qp.Organize.IC(n_g+n_l+i)>0 %all buildings have an initial air zone temperature state
        qp.beq(qp.Organize.IC(n_g+n_l+i)) = buildings(i).Tzone;
        qp.ub(qp.Organize.IC(n_g+n_l+i)) = buildings(i).Tzone+1;%+1 just to help solver find feasible
    end
end
for i = 1:1:n_fl
    if qp.Organize.IC(n_g+n_l+n_b+i)>0 %all fluid loops have an initial temperature state
        qp.beq(qp.Organize.IC(n_g+n_l+i)) = fluid_loop(i).fluid_temperature;
        qp.ub(qp.Organize.IC(n_g+n_l+i)) = fluid_loop(i).fluid_temperature+1;%+1 just to help solver find feasible
    end
end

qp = update_mat_demand(qp,gen,subnet,forecast,n_s);
qp = update_mat_costs(qp,gen,margin_cost,scale_cost,ec,dt);
qp = update_mat_build(qp,buildings,forecast,n_g,n_l,n_s);
if isfield(forecast,'Renewable')% subtract renewable generation 
    qp = update_mat_renew(qp,gen,forecast);
end
if isfield(subnet,'Hydro')
    qp = update_mat_hydro(qp,gen,subnet,forecast,date);
end 
if options.SpinReserve
    qp = update_mat_sr(qp,options,forecast,dt,n_g);
end
end%Ends function update_matrices

function qp = update_mat_demand(qp,gen,subnet,forecast,n_s)
%% update demands and storage self-discharge
n_g = length(gen);
network_names = fieldnames(subnet);
for nn = 1:1:length(network_names)
    net = network_names{nn};
    abrev = subnet.(net).abbreviation;
    qp.network.(abrev).nodes = length(subnet.(net).nodes);
    qp.constDemand.(abrev).req = zeros(qp.network.(abrev).nodes*n_s,1);
    qp.constDemand.(abrev).load = zeros(qp.network.(abrev).nodes*n_s,n_g);
    for n = 1:1:qp.network.(abrev).nodes %run through all the nodes in this network
        equip = subnet.(net).Equipment{n}; %equipment at this node
        req = qp.Organize.Balance.(net)(n);%balance at this node (t = 1)
        req = req:qp.Organize.t1Balances:((n_s-1)*qp.Organize.t1Balances+req);%balance at this node (t = 1:nS)
        load = nonzeros(qp.Organize.Demand.(net){n}); %loads at this node
        if ~isempty(load) %need this in case there is no field Forecast.Demand
            qp.beq(req) = sum(forecast.Demand.(abrev)(:,load),2); %multiple demands can be at the same node, or none
        end
        for j = 1:1:length(equip)
            k = equip(j);
            if ~isempty(strfind(gen(k).Type,'Storage'))
                if isfield(gen(k).QPform.output,abrev)
                    loss = (gen(k).QPform.Stor.SelfDischarge*gen(k).QPform.Stor.UsableSize*gen(k).QPform.Stor.DischEff);
                    qp.beq(req) = qp.beq(req) + loss; %account for self-discharge losses
                end
            end 
            if isfield(gen(k).QPform,'constDemand')  && isfield(gen(k).QPform.constDemand,abrev) && strcmp(qp.Organize.Fit,'B') %when using fit B record the constant electrical demands of the chillers
                qp.constDemand.(abrev).req((n-1)*n_s+1:n*n_s,1) = req;
                qp.constDemand.(abrev).load((n-1)*n_s+1:n*n_s,k) = gen(k).QPform.constDemand.(abrev);
            end
        end
    end
end  
end%Ends function update_mat_demand

function qp = update_mat_costs(qp,gen,margin_cost,scale_cost,ec,dt)
%update costs
n_g = length(gen);
n_s = length(dt);
H = diag(qp.H);
for i = 1:1:n_g
    if ~isempty(qp.Organize.States{i})
        states = qp.Organize.States{i};
        all_states = (0:qp.Organize.t1States:(n_s-1)*qp.Organize.t1States)'*ones(1,length(states))+ones(n_s,1)*states;
        switch gen(i).Type
            case 'Utility'
                qp.f(all_states(:,1)) = qp.f(all_states(:,1)).*scale_cost(:,i).*dt;
                if length(states)>1 && gen(i).VariableStruct.SellBackRate == -1
                    qp.f(all_states(:,2)) = qp.f(all_states(:,2)).*scale_cost(:,i).*dt; %sellback is a fixed percent of purchase costs (percent set when building matrices)
                else
                    %constant sellback rate was taken care of when building the matrices
                end
            case {'Electric Generator';'CHP Generator';'Chiller';'Heater';'Hydrogen Generator';'Cooling Tower';'Electrolyzer'}%all generators 
                for j = 1:1:length(states)
                    H(all_states(:,j)) = H(all_states(:,j)).*scale_cost(:,i).*dt; 
                    qp.f(all_states(:,j)) = qp.f(all_states(:,j)).*scale_cost(:,i).*dt; 
                end   
            case {'Thermal Storage'; 'Electric Storage';'Hydro Storage';'Hydrogen Storage';}
                if strcmp(gen(i).Type,'Hydro Storage')
                    s_end = (n_s-1)*qp.Organize.t1States + states(2); %SOC is second state
                else
                    s_end = (n_s-1)*qp.Organize.t1States + states(1);
                    if isfield(gen(i).QPform,'Y')%charging penalty state
                        qp.f(states(2):qp.Organize.t1States:s_end+1) = 1e-6; %small penalty on charging state so that if there are no other generators it keeps the charging constraint as an equality
                    end
                end

                stor_size = gen(i).QPform.Stor.UsableSize;
                if isfield(gen(i).QPform,'U') %has buffer
                    buff_size = gen(i).QPform.U.ub;
                else
                    buff_size = 0;
                end
                if isempty(ec) %full forecast optimization
                    S = fieldnames(gen(i).QPform.output);
                    max_value = 1.25*margin_cost.(S{1}).Max;
                    min_value = 0.5*margin_cost.(S{1}).Min;
                    a1 = -max_value; % fitting C = a1*SOC + 0.5*a2*SOC^2 so that dC/dSOC @ 0 = -max & dC/dSOC @ UB = -min
                    a2 = (max_value - min_value)/(stor_size);
                    H(s_end) = a2;%quadratic final value term loaded into SOC(t=nS)  %quadprog its solving C = 0.5*x'*H*x + f'*x
                    qp.f(s_end) = a1;%linear final value term loaded into SOC(t=nS)
                    if buff_size>0
                        qp.f(all_states(:,end-1)) = min_value;%this is the linear buffer term loaded into the lower & upper buffer
                        H(all_states(:,end-1)) = (2*max_value-min_value)/buff_size;%this is the quadratic buffer term loaded into the lower & upper buffer
                        qp.f(all_states(:,end)) = min_value;%this is the linear buffer term loaded into the lower & upper buffer
                        H(all_states(:,end)) = (2*max_value-min_value)/buff_size;%this is the quadratic buffer term loaded into the lower & upper buffer
                    end
                else %used in Online loop when there is a target EC determined by dispatch loop
                %final SOC deviation cost (quadratic cost for error = SOC(end) - EC
%                 if strcmp(Gen(i).Type,'Electric Storage')
%                     type = 'Electrical';
%                 elseif strcmp(Gen(i).Type,'Thermal Storage') && strcmp(Gen(i).Source,'Heat')
%                     type = 'DistrictHeat';
%                 elseif strcmp(Gen(i).Type,'Thermal Storage')
%                     type = 'DistrictCool';
%                 elseif strcmp(Gen(i).Type,'Hydro Storage')
%                     type = 'Hydro';
%                 end
%                 rows = QP.Organize.Inequalities{i};
%                 nr = length(states)/nS;
%                 PeakChargePower = Gen(i).QPform.Ramp.b(1);
%                 dSOC_10perc = .1*PeakChargePower*(Date(end)-Date(1)); %energy (kWh) if charging at 10%
%                 H(s_end) = -2*marginCost.(type).Min/dSOC_10perc;%quadratic final value term loaded into SOC(t=nS)  %factor of 2 because its solving C = 0.5*x'*H*x + f'*x
%                 QP.f(s_end) = -marginCost.(type).Min;%linear final value term loaded into SOC(t=nS)
%                 nIC = nnz(QP.Organize.IC(1:i)); %order in IC
%                 QP.beq(nIC) = IC(I)-EC(I);
%                 for t = 1:1:nS
%                     Xn = states((t-1)*nt+1);
%                     Rn = rows(t*nr);
%                     QP.lb(Xn) = -EC(i);%change lb so that SOC = 0 coresponds to EC
%                     QP.ub(Xn) = StorSize - EC(i);%change ub so that SOC = 0 coresponds to EC
%                     QP.b(Rn-1) = -BuffSize + EC(i);%change lb so that SOC = 0 coresponds to EC (adding EC because there is a -1 in front of SOC in this inequality)
%                     QP.b(Rn) = StorSize-BuffSize - EC(i);%change lb so that SOC = 0 coresponds to EC
%                 end       
                end    
        end
    end
end 
qp.H = diag(H);
end%End function update_mat_costs

function qp = update_mat_build(qp,buildings,forecast,n_g,n_l,n_s)
%Building inequalities
n_b = length(buildings);
for i = 1:1:n_b %all Forecasts are vectors of n_s steps
    states = qp.Organize.States{n_g+n_l+i};
    temperature_state = states(1):qp.Organize.t1States:(n_s-1)*qp.Organize.t1States + states(1);

    r = qp.Organize.Building.r(i):qp.Organize.t1ineq:(n_s-1)*qp.Organize.t1ineq + qp.Organize.Building.r(i);
    qp.b(r) = buildings(i).QPform.UA*forecast.Building.Tset_H(:,i);%heating inequality
    qp.b(r+1) = -buildings(i).QPform.UA*forecast.Building.Tset_C(:,i);%cooling inequality
    qp.b(r+2) = forecast.Building.Tmax(:,i);%upper buffer inequality
    qp.b(r+3) = -forecast.Building.Tmin(:,i);%lower buffer inequality
    
    req = qp.Organize.Building.Electrical.req(i):qp.Organize.t1Balances:(n_s-1)*qp.Organize.t1Balances + qp.Organize.Building.Electrical.req(i);
    qp.beq(req) = qp.beq(req) + forecast.Building.E0(:,i);%electric load
    %% update heating eqn
    req = qp.Organize.Building.DistrictHeat.req(i):qp.Organize.t1Balances:(n_s-1)*qp.Organize.t1Balances + qp.Organize.Building.DistrictHeat.req(i);
    qp.Organize.Building.H_Offset(:,i) = buildings(i).QPform.UA*(forecast.Building.Tzone(:,i)-forecast.Building.Tset_H(:,i)) - forecast.Building.H0(:,i);
    qp.beq(req) = qp.beq(req) - qp.Organize.Building.H_Offset(:,i);%Heating load
    qp.lb(temperature_state+1) = qp.lb(temperature_state+1) + qp.Organize.Building.H_Offset(:,i);%Ensures net building heating is greater than zero
    
    %% update cooling eqn
    req = qp.Organize.Building.DistrictCool.req(i):qp.Organize.t1Balances:(n_s-1)*qp.Organize.t1Balances + qp.Organize.Building.DistrictCool.req(i);
    qp.Organize.Building.C_Offset(:,i) = buildings(i).QPform.UA*(forecast.Building.Tset_C(:,i)-forecast.Building.Tzone(:,i)) - forecast.Building.C0(:,i);
    qp.beq(req) = qp.beq(req) - qp.Organize.Building.C_Offset(:,i);%Cooling load
    qp.lb(temperature_state+2) = qp.lb(temperature_state+2) + qp.Organize.Building.C_Offset(:,i);%Ensures net building cooling is greater than zero
end
end%Ends function update_mat_build

function qp = update_mat_renew(qp,gen,forecast)
qp.Renewable = forecast.Renewable;
n_s = length(forecast.Renewable(:,1));
n_g = length(gen);
for i = 1:1:n_g
    switch gen(i).Type
        case 'Solar'
            if isfield(gen(i).QPform,'Electrical')
                req = qp.Organize.Balance.Electrical(gen(i).QPform.Electrical.subnetNode);%balance at this node (t = 1)
            else
                req = qp.Organize.Balance.DirectCurrent(gen(i).QPform.DirectCurrent.subnetNode);%balance at this node (t = 1)
            end
            req = req:qp.Organize.t1Balances:((n_s-1)*qp.Organize.t1Balances+req);%balance at this node (t = 1:nS)
            qp.beq(req) = qp.beq(req) - forecast.Renewable(:,i); %put renewable generation into energy balance at correct node
    end
end
end%ends function update_mat_renew

function qp = update_mat_hydro(qp,gen,subnet,forecast,date)
%Updating lower and upper bound for Hydro resevoir SOC
qp.Organize.hydroSOCoffset = zeros(1,length(subnet.Hydro.nodes));
n_g = length(gen);
n_s = length(date)-1;
for n = 1:1:length(subnet.Hydro.nodes) %run through all the nodes in this network
    %% Node inflows (source/sink terms and upstream flow at time t-T ago if t<T)
    %need to be able to forecast sink/source
    req = qp.Organize.Balance.Hydro(n);%mass balance at this node (t = 1)
    req = req:qp.Organize.t1Balances:((n_s-1)*qp.Organize.t1Balances+req);%balance at this node (t = 1:nS)
    qp.beq(req)= qp.beq(req) - forecast.Hydro.InFlow(:,n); %source/sink term (also previous controlled releases that precede date(1))
end
for i = 1:1:n_g %Make sure bounds are between 0 and the maximum usable size for each generator
    switch gen(i).Type
        case 'Hydro Storage'
            n = gen(i).QPform.Hydro.subnetNode;%dam #
            if isfield(forecast,'WYForecast') && date(end)<forecast.wy_forecast{end}.Timestamp(end)
                qp.Organize.hydroSOCoffset(n) = interp1(forecast.wy_forecast{end}.Timestamp,forecast.wy_forecast{end}.hydroSOC(:,n),date(end));% re-center so zero is actual where this target is
                states = qp.Organize.States{i}; %States for Reservoirs
                soc_states = states(2):qp.Organize.t1States:(n_s-1)*qp.Organize.t1States + states(2); %All time steps for resevoir SOC
                range = max(0.1*gen(i).QPform.Stor.UsableSize,1+abs(gen(i).CurrentState(2)-qp.Organize.hydroSOCoffset(n)));
                qp.ub(soc_states) = min(range,gen(i).QPform.Stor.UsableSize-qp.Organize.hydroSOCoffset(n));
                qp.lb(soc_states) = max(-range,-qp.Organize.hydroSOCoffset(n));
            end
            soc_now = gen(i).CurrentState(2)-qp.Organize.hydroSOCoffset(n);
            qp.beq(qp.Organize.IC(i)+1) = soc_now;%SOC in reservior (2nd ic for dam)
            qp.ub(qp.Organize.IC(i)+1) = soc_now+1;%+1 just to help solver find feasible
            qp.lb(qp.Organize.IC(i)+1) = soc_now-1;%-1 just to help solver find feasible
    end
end
%Adding cost for Line Transfer Penalties to differentiate from spillway flow
if ~isempty(subnet.Electrical.lineNames)
    for k = 1:1:length(subnet.Electrical.lineNames)
        line = subnet.Electrical.lineNumber(k);
        states = qp.Organize.States{n_g+line};
        if length(states)>1 %bi-directional transfer with penalty states
            s2 = states(2):qp.Organize.t1States:((n_s-1)*qp.Organize.t1States+states(2));%penalty from a to b
            s3 = states(3):qp.Organize.t1States:((n_s-1)*qp.Organize.t1States+states(3));% penalty from b to a
            qp.f(s2) = 0.001; %Adding .1 cent/kWhr cost to every Electric line penalty to differentiate from spill flow
            qp.f(s3) = 0.001; %Adding .1 cent/kWhr cost to every Electric line penalty to differentiate from spill flow
        end
    end 
end
%Adding lower bounds for Hydro lines that have seasonality constraints
if ~isempty(subnet.Hydro.lineNames)
    for k = 1:1:length(subnet.Hydro.lineNames)
        line = subnet.Hydro.lineNumber(k);
        states = qp.Organize.States{n_g+line};
        t_states = states:qp.Organize.t1States:((n_s-1)*qp.Organize.t1States+states);
        months = datevec(date(2:end));
        months = months(:,2);
        qp.lb(t_states,1) = subnet.Hydro.lineMinimum(k);
%         qp.lb(t_states,1) = subnet.Hydro.lineMinimum{k}(months,1);
    end
end
end%Ends function update_mat_hydro


function qp = update_mat_sr(qp,options,forecast,dt,n_g)
n_s = length(dt);
H = diag(qp.H);

qp.b(qp.Organize.SpinReserve) = -forecast.SRtarget;% -shortfall + SRancillary - SR generators - SR storage <= -SR target
sr_short = qp.Organize.SpinReserveStates(:,n_g+1);%cumulative spinning reserve shortfall at t = 1 --> nS
all_sr_short = (0:qp.Organize.t1States:(n_s-1)*qp.Organize.t1States)'+ones(n_s,1)*sr_short;
%     sr_ancillary = QP.Organize.SpinReserveStates(:,nG+2);%Ancillary spinning reserve value (negative cost) at t = 1 --> nS
if options.SpinReservePerc>5 %more than 5% spinning reserve
    spin_cost = 2*dt./forecast.SRtarget;% -shortfall + SRancillary - SR generators - SR storage <= -SR target
else
    req = nonzeros(qp.Organize.Balance.Electrical);
    spin_cost = zeros(n_s,1);
    for t = 1:1:n_s
        spin_cost(t) = 2*dt(t)./(0.05*sum(qp.beq(req)));% -shortfall + SRancillary - SR generators - SR storage <= -SR target
        req = req + qp.Organize.t1Balances;
    end
end
H(all_sr_short) = spin_cost;%effectively $2 per kWh at point where shortfall = spin reserve percent*demand or $2 at 5%
qp.f(all_sr_short) = 0.05*dt; % $0.05 per kWh
qp.H = diag(H);
end%Ends function update_mat_sr