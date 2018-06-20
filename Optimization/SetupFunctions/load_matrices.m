function qp = load_matrices(gen,building,fluid_loop,subnet,optimoptions,op,dt)
%builds constant matrices for multi-time-step optimization
%Demands, initial conditions, and utility costs must updated prior to optimization
n_g = length(gen);
n_b = length(building);
n_fl = length(fluid_loop);
n_l = 0;
networkNames = fieldnames(subnet);
for net = 1:1:length(networkNames)
    n_l = n_l + length(subnet.(networkNames{net}).lineNames);
end
if ~isempty(dt)
    n_s = length(dt);
    organize.IC = zeros(n_g+n_l+n_b,1);
else
    n_s = 0;
end

organize.States=cell(1,n_g+n_l+n_b+n_fl);
organize.Dispatchable = false(1,n_g);
organize.Transmission = zeros(n_l,1);
organize.Hydro = [];
organize.Building.r = zeros(n_b,1);
organize.Building.req = zeros(n_b,1);
organize.fluid_loop = zeros(n_fl,1);
organize.dt = dt;
organize.Out_vs_State = cell(1,n_g+n_l+n_b+n_fl);
qp.organize = cell(n_s+1,n_g+n_l+n_b+n_fl);
qp.constCost = zeros(1,n_g);
qp.excessHeat = optimoptions.excessHeat;
qp.excessCool = optimoptions.excessCool;

if n_s>0
    [qp,organize,ic] = count_ic(gen,qp,organize,n_g,n_b,n_l,n_fl);
else
    ic = 0;
end
[qp,organize,lb,ub] = count_states(gen,subnet,optimoptions,qp,organize,n_g,n_b,n_l,n_fl,ic,op,n_s);

%% Expand by the number of time steps
% # of states per time step is xL 
% Order of states will repeat for nS time steps
organize.t1States = length(lb);
if n_s>0    
    total_states = n_s*organize.t1States+ic;
    for t= 2:n_s
        for n = 1:1:length(qp.organize(1,:))
            qp.organize{t+1,n} = qp.organize{t,n} + organize.t1States;
        end
    end
    H = zeros(total_states,1);
    qp.f = zeros(total_states,1);
    qp.lb = zeros(total_states,1);
    qp.ub = zeros(total_states,1);

    for t= 1:n_s
        qp.lb(ic+(t-1)*organize.t1States+1:ic+t*organize.t1States)  = lb;
        qp.ub(ic+(t-1)*organize.t1States+1:ic+t*organize.t1States)  = ub;
    end
else
    H = zeros(organize.t1States,1);
    qp.f = zeros(organize.t1States,1);
    qp.ub = ub;
    qp.lb = lb;
end

[organize,ec] = count_equalities(gen,subnet,optimoptions.SpinReserve,organize,ic,optimoptions.endSOC,n_s);
organize = count_inequalities(gen,subnet,optimoptions.SpinReserve,organize,n_b,n_s);
if n_s>0
    total_equalities = n_s*organize.t1Balances + ic + length(ec);
    qp.Aeq = spalloc(total_equalities,total_states,total_equalities*organize.t1States);
    qp.beq = zeros(total_equalities,1);
    qp.Aeq(1:ic,1:ic) = eye(ic); %initial condition identity matrix
    total_inequalities = n_s*organize.t1ineq;
    qp.A = spalloc(total_inequalities,total_states,total_inequalities*organize.t1States);
    qp.b = zeros(total_inequalities,1);
else
    qp.Aeq = zeros(organize.t1Balances,organize.t1States);
    qp.beq = zeros(organize.t1Balances,1);
    qp.A = zeros(organize.t1ineq,organize.t1States);
    qp.b = zeros(organize.t1ineq,1);
end

%% call sub functions to fill in these matrices now that they are sized
if n_s>0
    [qp.Aeq,qp.beq,H,qp.f] = generator_equalities(gen,subnet,qp.Aeq,qp.beq,H,qp.f,organize,dt,op);
    [qp.Aeq,qp.beq] = end_state_constraint(gen,optimoptions.endSOC, qp.Aeq,qp.beq,qp.organize,ec);
    [qp.A,qp.b] = generator_inequalities(gen,qp.A,qp.b,organize,dt);
    [qp.A,qp.b] = spin_reserve_inequalities(gen,subnet,optimoptions.SpinReserve,qp.A,qp.b,organize,dt);
else
    [qp,H] = equality_constraints_step(gen,subnet,qp,H,organize);
    qp.A = storage_inequalities_step(gen,qp.A,organize);
    qp.Aeq = hydro_equalities_step(gen,subnet,qp.Aeq,organize);
    [qp.Aeq,qp.beq,qp.A] = spin_reserve_constraints_step(gen,subnet,optimoptions.SpinReserve,qp.Aeq,qp.beq,qp.A,qp.ub,organize);
    total_states = organize.t1States;
end
qp.Aeq = line_equalities(subnet,qp.Aeq,organize,n_s,n_g);
qp.A = line_inequalities(subnet,qp.A,organize,n_s,n_g);
[qp.Aeq,qp.A,H,qp.f,organize] = building_constraints(building,qp.Aeq,qp.A,H,qp.f,organize,dt,n_g,n_l,n_s);
[qp.Aeq,qp.A,H,qp.f,organize] = fluid_loop_constraints(gen,fluid_loop,subnet,qp.Aeq,qp.A,H,qp.f,organize,dt,n_l,n_b,n_s);

%convert to sparse matrices
qp.H = spalloc(total_states,total_states,total_states);
for i = 1:1:total_states
    qp.H(i,i) = H(i);
end
qp.Aeq = sparse(qp.Aeq);
qp.A = sparse(qp.A);
qp.Organize = organize; %indices (rows and columns) associated with each generator, allowing generators to be removed later
end%Ends function build_matrices

function [qp,organize,ic] = count_ic(gen,qp,organize,n_g,n_b,n_l,n_ct)
% IC for each generator & storage
ic = 0; % row index of the initial condition constraint in the Aeq matrix and beq vector
for i = 1:1:n_g
    if isfield(gen(i).QPform,'Ramp') 
        ic = ic+1;%initial condition state
        qp.organize{1,i} = ic; %output state organized into matrix of time vs. generator (IC)   
        organize.IC(i) = ic;
    end
    if strcmp(gen(i).Type,'Hydro Storage')
        organize.Hydro(end+1) = i;
        ic = ic+1; %second IC for the SOC of the reservior (1st was for Power output)
    end
end
for i = 1:1:n_b %all buildings have an initial temperature state
    ic = ic+1;%initial condition state
    qp.organize{1,n_g+n_l+i} = ic; %output state organized into matrix of time vs. generator (IC)   
    organize.IC(n_g+n_l+i) = ic;
end
for i = 1:1:n_ct %all cooling tower water loops have an initial temperature state
    ic = ic+1;%initial condition state
    qp.organize{1,n_g+n_l+n_b+i} = ic; %output state organized into matrix of time vs. generator (IC)   
    organize.IC(n_g+n_l+n_b+i) = ic;
end
end%Ends function count_ic

function [qp,organize,lb,ub] = count_states(gen,subnet,optimoptions,qp,organize,n_g,n_b,n_l,n_ct,x_l,op,n_s)
%% organize the order of non-initial condition states 
% states for each generator/storage at t = 1
% spinning reserve for generator/storage at t = 1 if option is selected
% states for each transmission line/pipe at t = 1
% state for P_loss in heating & cooling energy balances
% state for spinning reserve shortfall (target - cumulative of all active generators) and SR provided as ancillary service
% repeat order of generators and lines for t = 2:nS
lb =[]; ub = [];
if optimoptions.SpinReserve
    organize.SpinReserveStates = zeros(1,n_g+2);
else
    organize.SpinReserveStates = [];
end
organize.Fit = op;
for i = 1:1:n_g
    gen_i = gen(i).QPform;
    if strcmp(op,'B') 
        [~,fit] = size(gen_i.states);% fit = 2 for generators with 2 different piecewise quadratics when Op = 'B'
    else
        fit = 1;
    end
    if ~isempty(gen_i.states)
        if n_s == 0 && isfield(gen_i,'Stor') %single step (convert storage to generator)
            if ismember(gen(i).Type,{'Hydro Storage'})
                s = 1;% Storage treated as generator with 1 state
                lb(end+1,1) = gen_i.X.lb; 
                ub(end+1,1) = gen_i.X.ub; % Upper bound is MaxGenCapacity
                if isfield(gen_i,'S') %spill flow
                    s = 2;
                    lb(end+1,1) = gen_i.S.lb; 
                    ub(end+1,1) = gen_i.S.ub; 
                end
            else
                if ~isfield(gen_i,'Y')
                    s = 1;% Ideal storage treated as generator with 1 states (addtl power output relative to 1st dispatch)
                    lb(end+1,1) = -gen_i.Ramp.b(1);
                    ub(end+1,1) = gen_i.Ramp.b(2);
                else
                    s = 2;% Storage treated as generator with 2 states (addtl power output relative to 1st dispatch & charging penalty)
                    lb(end+1:end+2,1) = [-gen_i.Ramp.b(1);0];
                    ub(end+1:end+2,1) = [gen_i.Ramp.b(2);(1/gen_i.Stor.ChargeEff-gen_i.Stor.DischEff)*gen_i.Ramp.b(2)];
                end
            end 
        else %dispatchable generator or utility
            states = gen_i.states(1:nnz(~cellfun('isempty',gen_i.states(:,fit))),fit);
            s = length(states);%generator with multiple states
            for k = 1:1:s
                lb(end+1,1) = gen_i.(states{k}).lb(fit);
                ub(end+1,1) = gen_i.(states{k}).ub(fit);
            end
            if ismember(gen(i).Type,{'Electric Generator';'CHP Generator';'Heater';'Chiller';'Cooling Tower';'Electrolyzer';'Hydrogen Generator';})
                if isfield(gen_i,'constCost') 
                     qp.constCost(i) = gen_i.constCost;
                end
                if isfield(gen(i).VariableStruct,'StartCost') && gen(i).VariableStruct.StartCost>0
                    organize.Dispatchable(i) = true;
                elseif isfield(gen_i,'constDemand')
                    organize.Dispatchable(i) = true;
                elseif qp.constCost(i)>0 || gen_i.(gen_i.states{1,end}).lb(end)>0 % first state has a non-zero lower bound for second optimization or a nonzero cost at an output of zero
                    organize.Dispatchable(i) = true;
                end
            end
        end
        qp.organize{min(n_s+1,2),i} = x_l+1:x_l+s; %output is sum of multiple states at each time step
        organize.States(i)= {x_l+1:x_l+s};
        if strcmp(gen(i).Type,'Utility') && strcmp(gen(i).Source,'Electricity') && s ==2
            organize.Out_vs_State{1,i} = [1,-1];
        elseif strcmp(gen(i).Type,'AC_DC')
            organize.Out_vs_State{1,i} = [1,-1];%Value is thus the total transfer from AC to DC
        elseif isfield(gen_i,'Stor')
            qp.organize{min(n_s+1,2),i} = x_l+1; %output state for storage is only SOC (Power in the case of Hydro)
            organize.Out_vs_State{1,i} = 1;
        else
            organize.Out_vs_State{1,i} = ones(1,s);
        end
        x_l = x_l + s;
        if optimoptions.SpinReserve
            include = {'Electric Generator';'CHP Generator';'Hydro Storage';'Electric Storage';'Hydrogen Generator';};
            if ismember(gen(i).Type,include)
                if strcmp(gen(i).Type,'Hydro Storage')
                    lb(end+1,1) = 0;
                    ub(end+1,1) = gen(i).VariableStruct.MaxGenCapacity;
                elseif strcmp(gen(i).Type,'Electric Storage')
                    lb(end+1,1) = -gen_i.Stor.PeakCharge;
                    ub(end+1,1) = gen_i.Stor.PeakDisch;
                else
                    lb(end+1,1) = 0;
                    ub(end+1,1) = gen(i).Size;
                end
                organize.SpinReserveStates(i) = x_l+1; %state of spinning reserve at time 1
                x_l = x_l + 1;
            end
        end
    end
end

[qp,organize,x_l,lb,ub] = count_trans_lines(subnet,qp,organize,x_l,lb,ub,n_g,n_s);
if optimoptions.SpinReserve && any(organize.SpinReserveStates(1:n_g)>0)
    organize.SpinReserveStates(n_g+1) = x_l+1; %spinning reserve shortfall at t=1
    organize.SpinReserveStates(n_g+2) = x_l+2; %SR provided as ancillary service at t=1
    x_l = x_l + 2; % Add states for reserve shortfall and ancillary service: shortfall state is equal to reserve target + reserve sold as ancillary service - actual spinning reserve
    lb(end+1:end+2,1) = 0;
    ub(end+1:end+2,1) = sum(ub(nonzeros(organize.SpinReserveStates(1:n_g))));
end

[organize.HeatVented,x_l,lb,ub] = vented_energy(gen,subnet,optimoptions.excessHeat,x_l,lb,ub,'H');
[organize.CoolVented,x_l,lb,ub] = vented_energy(gen,subnet,optimoptions.excessCool,x_l,lb,ub,'C');
for i = 1:1:n_b %all single zone buildings have 5 states
    organize.States(n_g+n_l+i) = {[x_l+1, x_l+2, x_l+3, x_l+4, x_l+5]};
    qp.organize{min(n_s+1,2),n_g+n_l+i} = x_l+1; %output state for building is temperature
    organize.Out_vs_State{1,n_g+n_l+i} = 1;
    x_l = x_l+5; % 5 states are: Temperature, Heating, Cooling, exceeding upper comfort zone, exceeding lower comfort zone
    lb(end+1:end+5,1) = [10; 0; 0; 0; 0;];%Dont let building below 10 degrees C
    ub(end+1:end+5,1) = [35; 1e5; 1e5; 10; 10;];%Building stays below 35 degrees C, and excess temperature is less than 10 degrees
end
for i = 1:1:n_ct %all cooling tower cold water loops have 1 state
    organize.States(n_g+n_l+n_b+i) = {x_l+1};
    qp.organize{min(n_s+1,2),n_g+n_l+n_b+i} = x_l+1; %output state for building is temperature
    organize.Out_vs_State{1,n_g+n_l+n_b+i} = 1;
    x_l = x_l+1; % 1 states: Temperature of return water
    lb(end+1,1) = 24;%Temperature in C
    ub(end+1,1) = 35;%Temperature in C
end
end%Ends function count_states

function [e_vented,xL,lb,ub] = vented_energy(gen,subnet,vent,xL,lb,ub,param)
%add states at each node to allow puroseful energy discharge (venting heat or cooling)
if strcmp(param,'H') 
    net = 'DistrictHeat';
elseif strcmp(param,'C')
    net = 'DistrictCool';
end
e_vented = [];
if any(strcmp(net,fieldnames(subnet))) && vent == 1
    %%find maximum heat production possible
    maxGen = 0;
    n_g = length(gen);
    for i = 1:1:n_g
        if isfield(gen(i).QPform.output,param)
            if ~isempty(strfind(gen(i).Type,'Storage'))
                maxGen = maxGen + gen(i).QPform.Stor.PeakDisch;
            else
                if strcmp(param,'H') && strcmp(gen(i).Type,'CHP Generator')
                    if isfield(gen(i).Output,'Electricity')
                        maxGen = maxGen + max(gen(i).Size*gen(i).Output.Capacity./gen(i).Output.Electricity.*gen(i).Output.Heat);
                    else
                        maxGen = maxGen + max(gen(i).Size*gen(i).Output.Capacity./gen(i).Output.DirectCurrent.*gen(i).Output.Heat);
                    end                    
                elseif gen(i).QPform.output.(param)(end)>0
                    maxGen = maxGen + max(max(gen(i).QPform.output.(param)))*gen(i).Size;
                end
            end
        end
    end
    %%assume heat can be lost any any node in the network that has a device producing heat
    n = length(subnet.(net).nodes);
    ventNodes = 0;
    e_vented =zeros(1,n); %matrix for the state associated with venting heat at each district heating node, at each time step
    for i = 1:1:n
        genI = subnet.(net).Equipment{i};%%identify generators at this node
        for j = 1:1:length(genI)
            if isfield(gen(genI(j)).QPform.output,'H')
                ventNodes = ventNodes + 1; % Add single state for heat that is ventd to make energy equality true
                e_vented(1,i) = (xL+ventNodes);
                break
            end
        end
    end
    xL = xL+ventNodes;
    lb(end+1:end+ventNodes,1) = 0;
    ub(end+1:end+ventNodes,1) = maxGen;
end
end%Ends function vented_energy

function [qp,organize,x_l,lb,ub] = count_trans_lines(subnet,qp,organize,x_l,lb,ub,n_g,n_s)
%states for transmission lines
network_names = fieldnames(subnet);
for net = 1:1:length(network_names)
    if ~isempty(subnet.(network_names{net}).lineNames)
        limit = subnet.(network_names{net}).lineLimit;
        if strcmp(network_names{net},'Hydro')
            eff = [];
            minimum = subnet.(network_names{net}).lineMinimum; 
        else
            eff = subnet.(network_names{net}).lineEff;
            minimum = zeros(length(limit(:,1)),1);
        end
        for i = 1:1:length(limit(:,1))
            line  = subnet.(network_names{net}).lineNumber(i);
            qp.organize{min(n_s+1,2),n_g+line} = x_l+1; %line state organized into matrix of time vs. line state
            organize.Out_vs_State{1,n_g+line} = 1;
            if isempty(eff) || length(eff(i,:))==1 || eff(i,2)==0 %uni-directional transfer, or river segment, 1 state for each line
                organize.States(n_g+line)= {x_l+1};
                lb(end+1,1) = minimum(i);
                ub(end+1,1) = limit(i,1);
                x_l = x_l + 1;
            else% bi-directional transfer, 3 states for each line (state of the line and penalty term in each direction)
                organize.States(n_g+line)= {[x_l+1, x_l+2, x_l+3]};
                lb(end+1:end+3,1) = [-limit(i,2),0,0];
                ub(end+1:end+3,1) = [limit(i,1), limit(i,1)*(1-eff(i,1)), limit(i,2)*(1-eff(i,2))];
                x_l = x_l + 3;
            end
        end
    end
end
end%Ends function count_trans_lines

function [organize,ec] = count_equalities(gen,subnet,sr,organize,ic,end_soc,n_s)
%% Equality Constraints
% IC for each generator and storage device
% Electric energy balance @ each Electric subNet node  at t = 1, including transmission lines
% Heat balance @ each DistricHeat subNet node at t = 1...
% Cooling balance @ each DistrictCool subNet node at t = 1
% Any generator link equalities (linking states within a generator)
% Repeat power balance equalities and link equalities at t = 2:nS
network_names = fieldnames(subnet);
req = ic; % row index of the Aeq matrix and beq vector

% The following puts together the energy balance equations
% 1 energy/mass balance equation for each subNet node
% Nodes were agregated if their line efficiencies were 1
for net = 1:1:length(network_names)
    n = length(subnet.(network_names{net}).nodes);
    organize.Balance.(network_names{net}) = zeros(n,1);
    organize.Demand.(network_names{net}) = cell(n,1);
    for i = 1:1:n
        req = req+1;%there is an energy/mass balance at this node
        organize.Balance.(network_names{net})(i) = req;
        %%note any demands at this node
        if ~isempty(subnet.(network_names{net}).Load{i})
            organize.Demand.(network_names{net})(i) = subnet.(network_names{net}).Load(i);
        end
    end
end
n_g = length(gen);
organize.Equalities = zeros(n_g,2);
for i = 1:1:n_g
    %link is a field if there is more than one state and the states are linked by an inequality or an equality
    if isfield(gen(i).QPform,'link') && isfield(gen(i).QPform.link,'eq')
        [m,~] = size(gen(i).QPform.link.eq);
        organize.Equalities(i,1:2) = [req+1, m]; 
        req = req+m;
    end
end
ec = [];
if n_s>0
    if ~strcmp(end_soc,'Flexible')
        for i = 1:1:n
            if isfield(gen(i).QPform,'Stor')
                ec(end+1) = i;%%Special case of final state-of-charge equal to initial state or fixed value
            end
        end
    end
else
    organize.SpinRow = zeros(1,n_g);
    if ismember('Hydro',network_names)
        organize.HydroInequalities = zeros(n,1);
    end
    include = {'Electric Generator';'CHP Generator';'Hydrogen Generator';'Electric Storage';'Hydro Storage';};
    for i = 1:1:n
        if strcmp(gen(i).Type,'Hydro Storage')
            organize.Hydro(end+1) = i;
            organize.HydroEqualities(i) = req+1;
            req = req+1;
        end
        if ismember('Electrical',network_names) && sr && ismember(gen(i).Type,include)%spinning reserve: no ramping constraint for single time step, so SR can be an equalities 
            organize.SpinRow(i) = req+1;%2 constraints i) ramping, 2) capacity
            req = req+1;
        end
    end
end
organize.t1Balances = req - ic;
end%Ends function equality_constraints

function organize = count_inequalities(gen,subnet,sr,organize,n_b,n_s)
%% Inequality Constraints
% 2 ramping constraints for each generator at each time step
% constraints for each Energy storage 
% 2 constraints for each bi-directional transmission line
% repeat all of the above for each time 2:nS
n_g = length(gen);
r = 0; % row index of the A matrix & b vector
network_names = fieldnames(subnet);
%agregate spinning reserve shortfall
%currently only implemented for electric power
if ismember('Electrical',fieldnames(subnet)) && sr
    organize.SpinReserve = r+1;
    r = r+1; %inequality for net spinning reserve (total of all idividual generator spinning reserves
    include = {'Electric Generator';'CHP Generator';'Electric Storage';'Hydro Storage';'Hydrogen Generator';};
end
if n_s>0
    %Ramping & Generator Inequalities
    organize.Ramping = zeros(1,n_g);
    organize.Inequalities = zeros(n_g,2);
    organize.SpinRow = zeros(1,n_g);
    for i = 1:1:n_g
        if isfield(gen(i).QPform,'Ramp') 
            organize.Ramping(i) = r+1; 
            r = r+2;
        end
        if isfield(gen(i).QPform,'link') && isfield(gen(i).QPform.link,'ineq')%link is a field if there is more than one state and the states are linked by an inequality or an equality
            [m,~] = size(gen(i).QPform.link.ineq);
            organize.Inequalities(i,1:2) = [r+1, m];
            r = r+m;
        end
        if ismember('Electrical',network_names) && sr && ismember(gen(i).Type,include)%spinning reserve inequalities 
            organize.SpinRow(i) = r+1;%2 constraints i) ramping, 2) capacity
            r = r+2;
        end
    end
else %single time step matrices
    %charging penalty eqn
    organize.Inequalities = zeros(1,n_g);
    for i = 1:1:n_g
        if isfield(gen(i).QPform,'Stor') && ~strcmp(gen(i).Type,'Hydro Storage')
            organize.Inequalities(i) = r+1;
            r = r+1;
        end
    end
end

%Transmission line inequalities (penalty terms)
for net = 1:1:length(network_names)
    if isfield(subnet.(network_names{net}),'lineEff') && ~isempty(subnet.(network_names{net}).lineEff)
        eff = subnet.(network_names{net}).lineEff;
        for i = 1:1:length(eff(:,1))
            if length(eff(i,:))==1 || eff(i,2)==0 %uni-directional transfer, 1 state for each line
                %do nothing, no inequalities linking penalty states
            else%bi-directional power transfer
                organize.Transmission(subnet.(network_names{net}).lineNumber(i)) = r+1; 
                r = r+2;
            end
        end
    end
end

for i = 1:1:n_b
    organize.Building.r(i) = r+1;
    r = r+4; % (4 * # of zones) inequalities per building
end

% number of generator inequality constraints & energy imbalances at each time step will be r
% there are 2 ramping constraints on each generator/storage
organize.t1ineq = r;
end%Ends function inequality_constraints

function [Aeq,beq,H,f] = generator_equalities(gen,subnet,Aeq,beq,H,f,organize,dt,op)
%this function loads generators/chillers, etc into the apropriate energy balance of Aeq
n_s = length(dt);
t1_balances = organize.t1Balances;
t1_states = organize.t1States;
network_names = fieldnames(subnet);
for net = 1:1:length(network_names)
    out = subnet.(network_names{net}).abbreviation;
    if strcmp(network_names{net},'Hydro')
        Aeq = hydro_equalities(gen,subnet,Aeq,organize,dt);
    else
        T = (1:1:n_s)';
        for n = 1:1:length(subnet.(network_names{net}).nodes)
            req = organize.Balance.(network_names{net})(n);
            equip = subnet.(network_names{net}).Equipment{n}; %only things in this list should have an output or input in this network (electric, heating, cooling
            for k = 1:1:length(equip)
                gen_i = gen(equip(k)).QPform;
                states = organize.States{equip(k)};
                %link is a field if there is more than one state and the states are linked by an inequality or an equality
                if isfield(gen_i,'link') && isfield(gen_i.link,'eq')
                    req2 = organize.Equalities(equip(k),1);
                    for j = 1:1:organize.Equalities(equip(k),2)
                        for t = 1:1:n_s
                            Aeq((t-1)*t1_balances+req2+(j-1),(t-1)*t1_states+states) = gen_i.link.eq(j,:);
                        end
                        beq((T-1)*t1_balances+req2+(j-1)) = gen_i.link.beq(j);
                    end
                end
                if isfield(gen_i,'Stor') && ~strcmp(gen(equip(k)).Type,'Hydro Storage')
                    eff = gen_i.Stor.DischEff;
                    for t = 1:1:n_s
                        Aeq((t-1)*t1_balances+req,(t-1)*t1_states+states(1)) = -eff/dt(t); %SOC at t
                        if ismember('Y',gen_i.states)
                            Aeq((t-1)*t1_balances+req,(t-1)*t1_states+states(2)) = -1/dt(t); %charging penalty at t
                        end
                        if t==1
                            Aeq(req,organize.IC(equip(k))) = eff/dt(t);%SOC at IC
                        else
                            Aeq((t-1)*t1_balances+req,(t-2)*t1_states+states(1)) = eff/dt(t); %SOC at t-1
                        end
                    end
                elseif ~isempty(states)
                    if strcmp(op,'B') 
                        [~,fit] = size(gen_i.states);% fit = 2 for generators with 2 different piecewise quadratics when Op = 'B'
                    else
                        fit = 1;
                    end
                    state_names = gen_i.states(1:nnz(~cellfun('isempty',gen_i.states(:,fit))),fit);
                    for j = 1:1:length(state_names)
                        H((T-1)*t1_states+states(j)) = gen_i.(state_names{j}).H(fit);%generator  cost fit
                        f((T-1)*t1_states+states(j)) = gen_i.(state_names{j}).f(fit);
                    end
                    if length(gen_i.output.(out)(1,:))>1
                        output = gen_i.output.(out)(1:length(state_names),fit);
                    else
                        output = gen_i.output.(out);
                    end
                    for t = 1:1:n_s
                        Aeq((t-1)*t1_balances+req,(t-1)*t1_states+states) = output';
                    end
                end
            end
            %%any heat loss term to balance equality
            if  ~isempty(organize.HeatVented) && strcmp('DistrictHeat',network_names{net}) && organize.HeatVented(n)~=0
                for t = 1:1:n_s
                    Aeq((t-1)*t1_balances+req,(t-1)*t1_states+organize.HeatVented(n)) = -1;
                end
            end
            %%any cool loss term to balance equality
            if  ~isempty(organize.CoolVented) && strcmp('DistrictCool',network_names{net}) && organize.CoolVented(n)~=0
                for t = 1:1:n_s
                    Aeq((t-1)*t1_balances+req,(t-1)*t1_states+organize.CoolVented(n)) = -1;
                end
            end
        end
    end
end
end%Ends function generator_equalities

function Aeq = hydro_equalities(gen,subnet,Aeq,organize,dt)
%This function loads the mass and energy balances of the hydro network
n_s = length(dt);
n_g = length(gen);
t1_balances = organize.t1Balances;
t1_states = organize.t1States;
downriver = {};
for i = 1:1:length(subnet.Hydro.lineNames)
    name = subnet.Hydro.lineNames{i};
    k = strfind(name,'_');
    downriver(end+1) = {name(k(2)+1:end)};% node names of the downriver node
end

for i = 1:1:n_g
    if strcmp(gen(i).Type,'Hydro Storage')
        states = organize.States{i};
        n = gen(i).QPform.Hydro.subnetNode;
        req = organize.Balance.Hydro(n);
        index_dr = gen(i).QPform.Hydro.down_river;%lines leaving this node, i.e. this node is the upriver node (should be at most 1)
        if any(strcmp(subnet.Hydro.nodes{n},downriver))
            index_ur = subnet.Hydro.lineNumber(subnet.Hydro.up_river{n});
        else
            index_ur = []; %no upstream nodes
        end
        for t = 1:1:n_s
            %% Mass balance
            Aeq((t-1)*t1_balances+req,(t-1)*t1_states+states(2)) = -12.1/dt(t); %SOC at t in acre-ft / hours converted to ft^3/s:  43560 ft^3/acre-ft * 1 hr / 3600s = 12.1 cfs per acre-ft/hr
            if t==1
                Aeq(req,organize.IC(i)+1) = 12.1/dt(t);%SOC at IC is 1 after dam power IC
            else
                Aeq((t-1)*t1_balances+req,(t-2)*t1_states+states(2)) = 12.1/dt(t);
            end
            %river segments flowing into this node
            for j = 1:1:length(index_ur)
                T = subnet.Hydro.lineTime(subnet.Hydro.up_river{n}(j));
                tt = sum(dt(1:t));
                if tt<=T                            
                    %Do nothing; the inflow rate will be updated in update_matrices
                elseif tt>T && tt<=T+dt(1)%between initial condition & first step
                    frac = (tt-T)/dt(1);%portion of flow from step 1, remainder from  SourceSink + (1-frac)*Inflow(t=1) : subtracted from beq in update matrices
                    Aeq((t-1)*t1_balances+req,organize.States{n_g+index_ur(j)})= frac; % Qupriver at step 1, 
                else
                    step = 2;
                    while tt>(T+sum(dt(1:step)))
                        step = step+1;
                    end
                    frac = (tt-(T+sum(dt(1:step-1))))/dt(step);%portion of flow from step, remainder from previous step
                    Aeq((t-1)*t1_balances+req,(step-1)*t1_states+organize.States{n_g+index_ur(j)})= frac; % Qupriver at t - floor(T)
                    Aeq((t-1)*t1_balances+req,(step-2)*t1_states+organize.States{n_g+index_ur(j)})= (1-frac); % Qupriver at t - ceil(T)
                end 
            end 
            %water flow out of the node
            Aeq((t-1)*t1_balances+req,(t-1)*t1_states+organize.States{n_g+index_dr})= -1; %Qdownriver,  
            %% convert power to outflow (other part of balance is done in Generator Equalities when doing the electric network, this is the - Qdownriver
            Aeq((t-1)*t1_balances+ organize.Equalities(i,1),(t-1)*t1_states+organize.States{n_g+index_dr})= -1; % Power/(eff*head*84.76) + Spill - Qdownriver = 0
        end 
    else
        %% add water district here
    end
end
end%Ends function hydro_equalities


function [Aeq,beq] = end_state_constraint(gen,end_soc,Aeq,beq,organize,ec)
n = length(beq) - length(ec);
for j = 1:1:length(ec)
    n = n+1;
    i = ec(j);
    if isnumeric(end_soc)%assume a vector of fixed end value conditions
        Aeq(n,organize{end,i}) = 1;
        beq(n,1) = end_soc(i);
    elseif strcmp(end_soc,'Initial')
        beq(n,1) = 0;
        if strcmp(gen(i).Type,'Hydro Storage')
            Aeq(n,organize{1,i}+1) = 1;
            Aeq(n,organize{end,i}+1) = -1;
        else%electric and thermal storage
            Aeq(n,organize{1,i}) = 1;
            Aeq(n,organize{end,i}) = -1;
        end
    end
end
end%Ends function end_state_constraint

function [A,b] = generator_inequalities(gen,A,b,organize,dt)
n_g = length(gen);
n_s = length(dt);
t1ineq = organize.t1ineq;
t1States = organize.t1States;
for i = 1:1:n_g
    gen_i = gen(i).QPform;
    states = organize.States{i};
    %Ramping 
    if organize.Ramping(i)>0
        %%if storage, ramping only affects 1st state
        if isfield(gen_i,'Stor')
            ramp_states = states(1);
        else
            ramp_states = states;
        end
        r = organize.Ramping(i);
        for t = 1:1:n_s
            A((t-1)*t1ineq+r,(t-1)*t1States+ramp_states) = 1/dt(t);%ramp up 
            A((t-1)*t1ineq+r+1,(t-1)*t1States+ramp_states) = -1/dt(t);%ramp down
            if t ==1 %ramping from initial condition
                A(r,organize.IC(i)) = -1/dt(t); %ramp up 
                A(r+1,organize.IC(i)) = 1/dt(t); %ramp down 
            else %condition at previous time step
                A((t-1)*t1ineq+r,(t-2)*t1States+ramp_states) = -1/dt(t); %ramp up 
                A((t-1)*t1ineq+r+1,(t-2)*t1States+ramp_states) = 1/dt(t); %ramp down 
            end
            b((t-1)*t1ineq+r,1) = gen_i.Ramp.b(1); %ramp up 
            b((t-1)*t1ineq+r+1,1) = gen_i.Ramp.b(2); %ramp down 
        end
    end
    %Inequalities constraints
    if organize.Inequalities(i,1)>0
        ineq_row = organize.Inequalities(i,1);
        for t = 1:1:n_s
            for k = 1:1:organize.Inequalities(i,2)
                if k == 1 && isfield(gen_i,'Stor') && ~strcmp(gen(i).Type,'Hydro Storage') && ismember('Y',gen_i.states)
                    A((t-1)*t1ineq+ineq_row,(t-1)*t1States+states(1)) = gen_i.link.ineq(1,1);  % SOC at t  
                    A((t-1)*t1ineq+ineq_row,(t-1)*t1States+states(2)) = gen_i.link.ineq(1,2);  % charging state at t (-1)
                    if t ==1 %SOC change from IC
                        A(ineq_row,organize.IC(i)) = -gen_i.link.ineq(1,1); % SOC at t-1
                    else
                        A((t-1)*t1ineq+ineq_row,(t-2)*t1States+states(1)) = -gen_i.link.ineq(1,1);  % SOC at t-1
                    end
                    loss = (gen_i.Stor.SelfDischarge*gen_i.Stor.UsableSize);
                    b((t-1)*t1ineq+ineq_row+(k-1)) = gen_i.link.bineq(k) - loss*(1-gen_i.Stor.ChargeEff);
                else
                    A((t-1)*t1ineq+ineq_row+(k-1),(t-1)*t1States+states) = gen_i.link.ineq(k,:); %typically buffer on SOC
                    b((t-1)*t1ineq+ineq_row+(k-1)) = gen_i.link.bineq(k);
                end
            end
        end
    end
end
end%Ends function generator_inequalities

function [qp,H] = equality_constraints_step(gen,subnet,qp,H,organize)
%this function loads generators/chillers, etc into the apropriate energy balance of Aeq
network_names = fieldnames(subnet);
for net = 1:1:length(network_names)
    if strcmp(network_names{net},'Hydro')
        %no Hydro mass balance
    else
        out = subnet.(network_names{net}).abbreviation;
        for n = 1:1:length(subnet.(network_names{net}).nodes)
            req = organize.Balance.(network_names{net})(n);
            equip = subnet.(network_names{net}).Equipment{n};
            for k = 1:1:length(equip)
                i = equip(k);
                gen_i = gen(i).QPform;
                states = organize.States{i};
                if isfield(gen(i).QPform,'Stor') 
                    if isfield(gen(i).QPform.output,out)    
                        qp.Aeq(req,states(1)) = gen(i).QPform.output.(out)(1);%additional power beyond 1st dispatch
                        if ~strcmp(gen(i).Type,'Hydro Storage') && isfield(gen(i).QPform,'Y') %has charging penalty
                            qp.Aeq(req,states(1)+1) = -1;%charging penalty
                        end
                    end
                else
                    %link is a field if there is more than one state and the states are linked by an inequality or an equality
                    if isfield(gen_i,'link') && isfield(gen_i.link,'eq')
                        req2 = organize.Equalities(i,1);
                        for j = 1:1:organize.Equalities(i,2)
                            qp.Aeq(req2+(j-1),states) = gen_i.link.eq(j,:);
                            qp.beq(req2+(j-1)) = gen_i.link.beq(j);
                        end
                    end
                    if ~isempty(gen_i.states)
                        [~,fit] = size(gen_i.states);% fit = 2 for generators with 2 different piecewise quadratics when Op = 'B'
                        state_names = gen_i.states(1:nnz(~cellfun('isempty',gen_i.states(:,fit))),fit);
                        n_states = length(state_names);
                        for j = 1:1:n_states
                            H(states(j)) = gen_i.(state_names{j}).H(fit);
                            qp.f(states(j)) = gen_i.(state_names{j}).f(fit);
                        end
                        if length(gen(i).QPform.output.(out)(1,:))>1
                            output = gen(i).QPform.output.(out)(1:n_states,fit)';
                        else
                            output = gen(i).QPform.output.(out);
                        end
                        qp.Aeq(req,states) = output;
                    end
                end
            end
            %%any heat loss term to balance equality
            if strcmp('DistrictHeat',network_names{net}) && qp.excessHeat == 1 && organize.HeatVented(n)~=0
                qp.Aeq(req,organize.HeatVented(n)) = -1;
            end
            if strcmp('DistrictCool',network_names{net}) && qp.excessCool == 1 && organize.CoolVented(n)~=0
                qp.Aeq(req,organize.CoolVented(n)) = -1;
            end
        end
    end
end
end%Ends function equality_constraints_step

function A = storage_inequalities_step(gen,A,organize)
n_g = length(gen);
for i = 1:1:n_g
    if isfield(gen(i).QPform,'Stor') && ~strcmp(gen(i).Type,'Hydro Storage')
        states = organize.States{i};
        r = organize.Inequalities(i);
        A(r,states(1)) = -(1/gen(i).QPform.Stor.ChargeEff-gen(i).QPform.Stor.DischEff);
        A(r,states(2)) = -1;
    end
end
end%ends function storage_inequalities_step

function Aeq = hydro_equalities_step(gen,subnet,Aeq,organize)
if isfield(subnet,'Hydro')
    n_g = length(gen);
    upriver = {};
    upLines = [];
    for i = 1:1:length(subnet.Hydro.lineNames)
        name = subnet.Hydro.lineNames{i};
        k = strfind(name,'_');
        upriver(end+1) = {name(1:k(1)-1)}; %node names of the upriver node (origin of line segment)
        upLines(end+1) = i;
    end
    %% equality:  AdditionalPower/(Hd*eff*84.67) + Spill Flow - Outflow = Nominal Power Flow(from 1st dispatch) 
        %This will allow imposition of instream flow constraints later
    for n = 1:1:length(subnet.Hydro.nodes)
        i = subnet.Hydro.Equipment{n};
        index_ur = upLines(strcmp(subnet.Hydro.nodes{n},upriver));%lines leaving this node, i.e. this node is the upriver node (should be at most 1)
        line_out = subnet.Hydro.lineNumber(index_ur);
        if strcmp(gen(i).Type,'Hydro Storage')
            req = organize.HydroEqualities(i);
            states = organize.States{i};
            Aeq(req,states(1)) = gen(i).QPform.Stor.Power2Flow;
            Aeq(req,organize.States{n_g+line_out})= -1; %Qdownriver (in excess of original solution from multi-time step),  
            if isfield(gen(i).QPform,'S') %spill flow
                Aeq(req,states(2)) = 1;
            end
        end
    end
end
end%Ends function hydro_equalities_step

function Aeq = line_equalities(subnet,Aeq,organize,n_s,n_g)
%this function loads generators/chillers, etc into the apropriate energy balance of Aeq
if n_s == 0
    n_s = 1;
end
t1_balances = organize.t1Balances;
t1_states = organize.t1States;
network_names = fieldnames(subnet);
for net = 1:1:length(network_names)
    if ~strcmp(network_names{net},'Hydro') %hydro is done later because of the time of transfer of the lines
        node_names = cell(length(subnet.(network_names{net}).nodes),1);
        for n = 1:1:length(node_names)
            agregated_nodes = subnet.(network_names{net}).nodes{n};
            node_names(n) = agregated_nodes(1);
        end
        for k = 1:1:length(subnet.(network_names{net}).lineNames)
            line = subnet.(network_names{net}).lineNames{k};
            linenum = subnet.(network_names{net}).lineNumber(k);
            r = strfind(line,'_');
            origin = line(1:r(1)-1);
            destination = line(r(2)+1:end);
            n1 = find(strcmp(origin,node_names),1,'first');
            linestates = organize.States{n_g+linenum};
            req = organize.Balance.(network_names{net})(n1);
            if length(linestates)==1 %uni-directional transfer
                for t = 1:1:n_s
                    Aeq((t-1)*t1_balances+req,(t-1)*t1_states+linestates) = -1;
                end
            else
                for t = 1:1:n_s
                    Aeq((t-1)*t1_balances+req,(t-1)*t1_states+linestates) = [-1,0,-1]; %sending node (this node, i)-->(connected node), is positive, thus positive transmission is power leaving the node, the penalty from b->a is power not seen at a
                end
            end
            n2 = find(strcmp(destination,node_names),1,'first');
            req = organize.Balance.(network_names{net})(n2);
            if length(linestates)==1 %uni-directional transfer
                for t = 1:1:n_s
                    Aeq((t-1)*t1_balances+req,(t-1)*t1_states+linestates) = subnet.(network_names{net}).lineEff(k);
                end
            else
                for t = 1:1:n_s
                    Aeq((t-1)*t1_balances+req,(t-1)*t1_states+linestates) = [1,-1,0];%receiving node (connected node)-->(this node, i), is positive, thus positive power is power entering the node, the penalty from a->b is power not seen at b
                end
            end
        end
    end
end
end%Ends function line_equalities

function A = line_inequalities(subnet,A,organize,n_s,n_g)
%Transmission
%%no connection to previous or later time steps, and no dependence on step size. 
if n_s == 0
    n_s = 1;
end
t1_ineq = organize.t1ineq;
t1_states = organize.t1States;
network_names = fieldnames(subnet);
for net = 1:1:length(network_names)
    if isfield(subnet.(network_names{net}),'lineEff')&& ~isempty(subnet.(network_names{net}).lineEff)
        eff = subnet.(network_names{net}).lineEff;
        for i = 1:1:length(eff(:,1))
            line = subnet.(network_names{net}).lineNumber(i);
            line_row = organize.Transmission(line);
            states = organize.States{n_g+line};
            if~isempty(line_row)
                for t = 1:1:n_s
                    A((t-1)*t1_ineq+line_row,(t-1)*t1_states+states) = [(1-eff(i,1)), -1, 0];% Pab*(1-efficiency) < penalty a to b
                    A((t-1)*t1_ineq+line_row+1,(t-1)*t1_states+states) = [-(1-eff(i,2)), 0, -1];% -Pab*(1-efficiency) < penalty b to a
                end
            end
        end
    end
end
end%Ends function line_inequalities

function [A,b] = spin_reserve_inequalities(gen,subnet,sr,A,b,organize,dt)
%spinning reserve inequalities (sum all spinning reserves & individual spinning reserves) 
%currently only implemented for electric power
if ismember('Electrical',fieldnames(subnet)) && sr
    n_g = length(gen);
    n_s = length(dt);
    t1_ineq = organize.t1ineq;
    t1_states = organize.t1States;
    %total spinning reserve shortfall at each time step
    sr_states = nonzeros(organize.SpinReserveStates(1:n_g+1));
    sr_ancillary = organize.SpinReserveStates(n_g+2);
    for t = 1:1:n_s
        A((t-1)*t1_ineq+organize.SpinReserve,(t-1)*t1_states+sr_states) = -1; %Inequality for spinning reserve shortfall:  -(shortfall) - sum(SR(i)) + SR ancillary <= -SR target
        A((t-1)*t1_ineq+organize.SpinReserve,(t-1)*t1_states+sr_ancillary) = 1; %Inequality for spinning reserve shortfall:  -(shortfall) - sum(SR(i)) + SR ancillary <= -SR target
    end
    %individual spinning reserve at each time step (limited by ramp rate and peak gen capacity)
    include = {'Electric Generator';'CHP Generator';'Hydro Storage';'Hydrogen Generator'};
    for i = 1:1:n_g
        if ismember(gen(i).Type,include)
            sr_state = organize.SpinReserveStates(i);
            states = organize.States{i};
            if strcmp(gen(i).Type,'Hydro Storage')
                states = states(1);%power output only associated with 1st state
            end
            r = organize.SpinRow(i);
            for t = 1:1:n_s
                A((t-1)*t1_ineq+r,(t-1)*t1_states+sr_state) = 1/dt(t); %SR + power(t) - power(t-1)<= ramprate*dt
                A((t-1)*t1_ineq+r,(t-1)*t1_states+states) = 1/dt(t);
                if t ==1 %ramping from IC
                    A(r,organize.IC(i)) = -1/dt(t); % Power at t-1
                else
                    A((t-1)*t1_ineq+r,(t-2)*t1_states+states) = -1/dt(t); % Power at t-1
                end
                b((t-1)*t1_ineq+r) = gen(i).QPform.Ramp.b(1);%ramp up constraint

                A((t-1)*t1_ineq+r+1,(t-1)*t1_states+sr_state) = 1; %SR + power <= Size
                A((t-1)*t1_ineq+r+1,(t-1)*t1_states+states) = 1;
                b((t-1)*t1_ineq+r+1) = gen(i).Size; %max capacity constraint
            end
        end
    end
    %individual spinning reserve at each time step (limited by discharge rate and storage capacity)
    include = {'Electric Storage';};
    for i = 1:1:n_g
        if ismember(gen(i).Type,include)
            sr_state = organize.SpinReserveStates(i);
            states = organize.States{i};
            soc = states(1);
            r = organize.SpinRow{i};
            for t = 1:1:n_s
                eff = gen(i).QPform.Stor.DischEff;
                A((t-1)*t1_ineq+r,(t-1)*t1_states+sr_state) = 1; %SR + eff*(SOC(t-1) - SOC(t))/dt <= peak discharge
                A((t-1)*t1_ineq+r,(t-1)*t1_states+soc) = -eff/dt(t);
                A((t-1)*t1_ineq+r+1,(t-1)*t1_states+sr_state) = 1; %SR - SOC(t-1)/dt <= 0
                if t ==1 %SOC change from IC
                    A(r,organize.IC(i)) = eff/dt(t); % SOC at t-1
                    A(r+1,organize.IC(i)) = -eff/dt(t); % SOC at t-1
                else
                    A((t-1)*t1_ineq+r,(t-2)*t1_states+soc) = eff/dt(t); % SOC at t-1
                    A((t-1)*t1_ineq+r+1,(t-1)*t1_states+soc) = -eff/dt(t); % SOC at t-1
                end
                b((t-1)*t1_ineq+r) = gen(i).QPform.Ramp.b(2);%peak discharge constraint
            end
        end
    end
end
end%Ends function spin_reserve_inequalities

function [Aeq,beq,A] = spin_reserve_constraints_step(gen,subnet,sr,Aeq,beq,A,ub,organize)
%spinning reserve equalities: individual spinning reserves
%currently only implemented for electric power
if ismember('Electrical',fieldnames(subnet)) && sr
    n_g = length(gen);
    %individual spinning reserve at each time step (limited by ramp rate and peak gen capacity)
    include = {'Electric Generator';'CHP Generator';'Hydro Storage';'Electric Storage';'Hydrogen Generator';};
    for i = 1:1:n_g
        if ismember(gen(i).Type,include)
            states = organize.States{i};
            if strcmp(gen(i).Type,'Hydro Storage')
                states = states(1);%power output only associated with 1st state
            end
            req = organize.SpinRow(i);
            Aeq(req+1,states) = 1;
            beq(req+1) = sum(ub(states(1:end-1)));%max power
        end
    end
    A(organize.SpinReserve,nonzeros(organize.SpinReserveStates(1,1:n_g+1))) = -1; %Inequality for spinning reserve shortfall:  -(shortfall) - sum(SR(i)) + SR ancillary <= -SR target
    A(organize.SpinReserve,nonzeros(organize.SpinReserveStates(1,n_g+2))) = 1; %Inequality for spinning reserve shortfall:  -(shortfall) - sum(SR(i)) + SR ancillary <= -SR target
end
end%Ends function spin_reserve_constraints_step

function [Aeq,A,H,f,organize] = building_constraints(building,Aeq,A,H,f,organize,dt,n_g,n_l,n_s)
%States representing building are: Tzone, Hating, Cooling, Excess T, shortfall T
n_b = length(building);
if n_s>0
    t1_balances = organize.t1Balances;
    t1_ineq = organize.t1ineq;
    t1_states = organize.t1States;
    for i = 1:1:n_b
        states = organize.States{n_g+n_l+i};
        %Electric equality
        req = organize.Balance.Electrical(building(i).QPform.Electrical.subnetNode);
        organize.Building.Electrical.req(i) = req;
        for t = 1:1:n_s
            Aeq((t-1)*t1_balances+req,(t-1)*t1_states+states(2)) = -building(i).QPform.H2E;% Electricity = E0 + H2E*(Heating-H0) + C2E*(Cooling-C0)
            Aeq((t-1)*t1_balances+req,(t-1)*t1_states+states(3)) = -building(i).QPform.C2E;% Electricity = E0 + H2E*(Heating-H0) + C2E*(Cooling-C0)
        end
        %Heating equality
        req = organize.Balance.DistrictHeat(building(i).QPform.DistrictHeat.subnetNode);
        organize.Building.DistrictHeat.req(i) = req;
        for t = 1:1:n_s
            Aeq((t-1)*t1_balances+req,(t-1)*t1_states+states(2)) = -1;%subtract building heating needs from heating energy balance
        end
        %Cooling equality
        req = organize.Balance.DistrictCool(building(i).QPform.DistrictCool.subnetNode);
        organize.Building.DistrictCool.req(i) = req;
        for t = 1:1:n_s
            Aeq((t-1)*t1_balances+req,(t-1)*t1_states+states(3)) = -1;%subtract building cooling needs from cooling energy balance
        end
        r = organize.Building.r(i);

        for t = 1:1:n_s
            %heating inequality % H = H0 + Hbar where H_bar>=  UA*(Ti-Tset_H) + Cap*(Ti - T(i-1))/dt where dt is in seconds
            A((t-1)*t1_ineq+r,(t-1)*t1_states+states(1)) = (building(i).QPform.UA+building(i).QPform.Cap/(3600*dt(t)));
            A((t-1)*t1_ineq+r,(t-1)*t1_states+states(2)) = -1;
            if t == 1
                A(r,organize.IC(n_g+n_l+i)) = -building(i).QPform.Cap/(3600*dt(t));
            else
                A((t-1)*t1_ineq+r,(t-2)*t1_states+states(1)) = -building(i).QPform.Cap/(3600*dt(t));
            end
            %Cooling inequality % Cooling = C0 + C_bar where C_bar>= UA*(Tset_C-Ti) + Cap*(T(i-1)-Ti)/dt where dt is in seconds
            A((t-1)*t1_ineq+r+1,(t-1)*t1_states+states(1)) = -building(i).QPform.UA - building(i).QPform.Cap/(3600*dt(t));
            A((t-1)*t1_ineq+r+1,(t-1)*t1_states+states(3)) = -1;
            if t == 1
                A(r+1,organize.IC(n_g+n_l+i)) = -building(i).QPform.Cap/(3600*dt(t));
            else
                A((t-1)*t1_ineq+r+1,(t-2)*t1_states+states(1)) = -building(i).QPform.Cap/(3600*dt(t));% Cooling>= C0 + UA*(Tset_C-Ti) + Cap*(T(i-1)-Ti)/dt where dt is in seconds
            end
            %Upper bound inequality
            A((t-1)*t1_ineq+r+2,(t-1)*t1_states+states(1)) = 1;% Upper buffer >= Ti - (Tset + Comfort width/2)
            A((t-1)*t1_ineq+r+2,(t-1)*t1_states+states(4)) = -1;% Upper buffer >= Ti - (Tset + Comfort width/2)
            %Lower bound inequality
            A((t-1)*t1_ineq+r+3,(t-1)*t1_states+states(1)) = -1;% Lower buffer >= (Tset - Comfort width/2) - Ti
            A((t-1)*t1_ineq+r+3,(t-1)*t1_states+states(5)) = -1;% Lower buffer >= (Tset - Comfort width/2) - Ti
            %Cost penalty for exceeding temperature bounds
            H((t-1)*t1_states+states(4)) = dt(t)*building(i).QPform.Discomfort;
            H((t-1)*t1_states+states(5)) = dt(t)*building(i).QPform.Discomfort;
            f((t-1)*t1_states+states(4)) = dt(t)*building(i).QPform.Discomfort;
            f((t-1)*t1_states+states(5)) = dt(t)*building(i).QPform.Discomfort;
        end
    end
else
    for i = 1:1:n_b
        states = organize.States{n_g+n_l+i};
        %Electric equality %put demands into electrical (heating and cooling) balances
        req = organize.Balance.Electrical(building(i).QPform.Electrical.subnetNode);
        organize.Building.Electrical.req(i) = req;
        Aeq(req,states(2)) = -building(i).QPform.H2E;% Electricity = E0 + H2E*(Heating-H0) + C2E*(Cooling-C0)
        Aeq(req,states(3)) = -building(i).QPform.C2E;% Electricity = E0 + H2E*(Heating-H0) + C2E*(Cooling-C0)
        %Heating equality
        req = organize.Balance.DistrictHeat(building(i).QPform.DistrictHeat.subnetNode);
        organize.Building.DistrictHeat.req(i) = req;
        Aeq(req,states(2)) = -1;%subtract building heating needs from heating energy balance
        %Cooling equality
        req = organize.Balance.DistrictCool(building(i).QPform.DistrictCool.subnetNode);
        organize.Building.DistrictCool.req(i) = req;
        Aeq(req,states(3)) = -1;%subtract building cooling needs from cooling energy balance
        r = organize.Building.r(i);
        %heating inequality % H = H0 + Hbar where H_bar>=  UA*(Ti-Tset_H) + Cap*(Ti - T(i-1))/dt where dt is in seconds
        % Done in update1Step because of dt: % A(r,states(1)) = (Building(i).QPform.UA+Building(i).QPform.Cap/dt);
        A(r,states(2)) = -1;
        %Cooling inequality % Cooling = C0 + C_bar where C_bar>= UA*(Tset_C-Ti) + Cap*(T(i-1)-Ti)/dt where dt is in seconds
        % Done in update1Step because of dt: % A(r+1,states(1)) = (Building(i).QPform.UA-Building(i).QPform.Cap/dt);
        A(r+1,states(3)) = -1;
        %Upper bound inequality
        A(r+2,states(1)) = 1;% Upper buffer >= Ti - (Tset + Comfort width/2)
        A(r+2,states(4)) = -1;% Upper buffer >= Ti - (Tset + Comfort width/2)
        %Lower bound inequality
        A(r+3,states(1)) = -1;% Lower buffer >= (Tset - Comfort width/2) - Ti
        A(r+3,states(5)) = -1;% Lower buffer >= (Tset - Comfort width/2) - Ti
        %Cost penalty for exceeding temperature bounds
        H(states(4)) = building(i).QPform.Discomfort;
        H(states(5)) = building(i).QPform.Discomfort;
        f(states(4)) = building(i).QPform.Discomfort;
        f(states(5)) = building(i).QPform.Discomfort;
    end
end
end%Ends function building_constraints

function [Aeq,A,H,f,organize] = fluid_loop_constraints(gen,fluid_loop,subnet,Aeq,A,H,f,organize,dt,n_l,n_b,n_s)
n_g = length(gen);
n_ct = length(fluid_loop);
if n_s>0
    t1Balances = organize.t1Balances;
    t1States = organize.t1States;
    for i = 1:1:n_ct
        %1st order temperature model: 0 = -T(k) + T(k-1) + dt/Capacitance*(energy balance)
        state = organize.States{n_g+n_l+n_b+i};
        req = organize.Balance.CoolingWater(i);
        organize.fluid_loop(i) = req;
        capacitance = fluid_loop(i).fluid_capacity*fluid_loop(i).fluid_capacitance; %Water capacity in kg and thermal capacitance in kJ/kg*K to get kJ/K
        for t = 1:1:n_s%energy balance for chillers & cooling tower fans already put in this row of Aeq
            Aeq((t-1)*t1Balances+req,(t-1)*t1States+state) = -capacitance/(3600*dt(t)); %-T(k)
            if t == 1
                Aeq(req,organize.IC(n_g+n_l+n_b+i)) = capacitance/(3600*dt(t));%T(k-1)
            else
                Aeq((t-1)*t1Balances+req,(t-2)*t1States+state) = capacitance/(3600*dt(t));%T(k-1)
            end
        end        
        delta_temperature = fluid_loop(i).nominal_supply_temperature - fluid_loop(i).nominal_return_temperature;
        pump_power = fluid_loop(i).pump_power_per_kgs/(delta_temperature*fluid_loop(i).fluid_capacitance);
        equip = subnet.CoolingWater.Equipment{i};
        for k = 1:1:length(equip)
            j = equip(k);
            if strcmp(gen(j).Type,'Chiller')
                states = organize.States{j};
                req = organize.Balance.Electrical(gen(i).QPform.Electrical.subnetNode);%Electric equality
                for t = 1:1:n_s
                    Aeq((t-1)*t1Balances+req,(t-1)*t1States+states) = Aeq((t-1)*t1Balances+req,(t-1)*t1States+states) - pump_power;%Add pump loads for condensor/cooling tower loop
                end
            end
        end
    end
else
    %single step optimization
    for i = 1:1:n_ct
        %1st order temperature model: T(k) = T(k-1) + dt/Capacitance*(energy balance)
        % done in update1Step because of dt: % Aeq(req,chiller/ fan states) = (HeatRejected)/Building(i).QPform.Cap*dt;
        req = organize.Balance.CoolingWater(i);
        organize.fluid_loop(i) = req;
        delta_temperature = fluid_loop(i).nominal_supply_temperature - fluid_loop(i).nominal_return_temperature;
        pump_power = fluid_loop(i).pump_power_per_kgs/(delta_temperature*fluid_loop(i).fluid_capacitance);
        equip = subnet.CoolingWater.Equipment{i};
        for k = 1:1:length(equip)
            j = equip(k);
            if strcmp(gen(j).Type,'Chiller')
                states = organize.States{j};
                req = organize.Balance.Electrical(gen(i).QPform.Electrical.subnetNode);%Electric equality
                Aeq(req,states) = Aeq(req,states) - pump_power;%Add pump loads for condensor/cooling tower loop
            end
        end
    end
end
end%Ends function fluid_loop_constraints