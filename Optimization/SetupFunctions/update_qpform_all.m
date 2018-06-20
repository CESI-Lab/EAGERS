function [gen, buildings] = update_qpform_all(gen,buildings,network,scaletime)% updates the QPform field in all generators and buildings
%% this function identifies the values that will be used to represent each generator in the quadratic optimizations
%% Network names, their corresponding abbreviation, what they represent
% Electrical --- E  --- standard 480V AC electrical network 
% DistrictHeat --- H --- standard 80C supply heating bus
% DistrictCool --- C --- standard 4C supply cooling bus
% Hydro        --- W --- River network with reservoirs and dams
% DirectCurrent --- DC --- 48V DC electrical network
% CoolingWater --- CW --- Water circulated between chillers and cooling towers
% Transmission1 --- E1 --- 230kV electric transmission (E2, E3, etc can be additional voltage levels
% Hydrogen     --- Hy --- Gaseous hydrogen stream
% LiqHydrogen  --- LH2 --- Liquid hydrogen
% Heating2     --- H2 --- Heat a different temperature than DistrictHeat (H3, H4... as needed)
%%-----%%%
n_g = length(gen);
n = 5; % # of segments in piecewise quadratic fits
for i = 1:1:n_g
    switch gen(i).Type
        case 'Utility'
            gen(i).VariableStruct.dX_dt = inf;
            gen(i).QPform = load_utility(gen(i).VariableStruct,gen(i).Source);
        case {'Electric Generator';'CHP Generator';'Chiller';'Heater';'Cooling Tower';'Electrolyzer';'Hydrogen Generator'}
            gen(i).QPform = load_piecewise(gen(i),n);
            gen(i).VariableStruct.dX_dt = ss_response(gen(i),[])/scaletime;
            gen(i).QPform.Ramp.b = gen(i).VariableStruct.dX_dt*[1;1]; %-output1+output2=ramp up %output1-output2=-rampdown
        case 'Solar'
            gen(i).QPform.states = [];
            if ~isfield(network,'DirectCurrent')
                gen(i).QPform.output.E = 1;%there are no states or outputs for solar because renewable outputs are handled on the demand side
            else 
                gen(i).QPform.output.DC = 1;
            end
        case 'Wind'
            
        case {'Electric Storage'}
            [gen(i).QPform,gen(i).VariableStruct.dX_dt] = load_storage(gen(i),scaletime);
            if ~isfield(network,'DirectCurrent')
                gen(i).QPform.output = [];
                gen(i).QPform.output.E = 1; 
            end
        case {'Thermal Storage';'Hydrogen Storage'}
            [gen(i).QPform,gen(i).VariableStruct.dX_dt] = load_storage(gen(i),scaletime);
        case 'Hydro Storage'
            gen(i).QPform = load_hydro_storage(gen(i),scaletime);
            gen(i).VariableStruct.dX_dt = gen(i).VariableStruct.RampUp/scaletime;
        case 'AC_DC'
            gen(i).QPform = load_ac_dc(gen(i));
    end
    gen(i).CurrentState = [];
    gen(i).Status = true;
end 
%% this loads a representation of a building a single zone temperature 
% Within the optimization 5 states are used:
% #1 represents the air zone temperature setpoint
% #2 represents the heating required to achieve that zone temperature setpoint
% #3 represents the cooling required to achieve that zone temperature setpoint
% #4 represents the temperature in excess of a comfort range (this is penalized)
% #5 represents the temperature below a comfort range (this is penalized)
n_b = length(buildings);
for i = 1:1:n_b
    f = .1;%fraction of wall R-value attributed to between wall and zone, remainder of R-value is between ambient and wall (needs to be 0.1 for small office)
    ua_window = buildings(i).VariableStruct.WallArea*buildings(i).VariableStruct.WindowWallRatio*buildings(i).VariableStruct.WindowUvalue/1000;
    ua_wall = (buildings(i).VariableStruct.WallArea*(1-buildings(i).VariableStruct.WindowWallRatio)/buildings(i).VariableStruct.WallRvalue + buildings(i).VariableStruct.RoofArea/buildings(i).VariableStruct.RoofRvalue)/1000;% Transmittance in kW/K, see the following about whole wall R-Value http://web.ornl.gov/sci/buildings/docs/Thermal-Performance-and-Wall-Ratings.pdf
    air_capacitance = 2*buildings(i).Area;

    qp_form.states = {'T';'H';'C';'U';'L'};%Temperature, Heating, Cooling, exceeding upper comfort zone, exceeding lower comfort zone
    qp_form.UA = ua_window+1/f*ua_wall;
    qp_form.Cap = air_capacitance;
    qp_form.H2E = 1/buildings(i).VariableStruct.COP_H;
    qp_form.C2E = 1/buildings(i).VariableStruct.COP_C;
    qp_form.Cooling = false; %initial condition is that heating and cooling are converted to electricity, this can change in locate_buildings if there are local chillers or heaters
    qp_form.Heating = false;
    qp_form.Discomfort = buildings(i).Area*buildings(i).VariableStruct.occupancy; %cost of exceeding comfort band
    buildings(i).QPform = qp_form;
    buildings(i).Tzone = 20;
    buildings(i).Twall = 20;
    buildings(i).Timestamp = 0;
end 
end%Ends function load_generator

function qp_form = load_ac_dc(gen)
qp_form.states = {'A';'B'};%first state is AC power transfered to DC power, second state is DC power to AC power
if gen.VariableStruct.DC_to_AC_eff == 1 && gen.VariableStruct.AC_to_DC_eff==1
    gen.VariableStruct.AC_to_DC_eff = 1-1e-4;%prevent infinite back/forward transfer in ideal AC_DC
end
qp_form.output.E = [-1;gen.VariableStruct.DC_to_AC_eff];
qp_form.output.DC = [gen.VariableStruct.AC_to_DC_eff;-1];
qp_form.A.ub = gen.VariableStruct.Capacity;
qp_form.A.lb = 0;
qp_form.B.ub = gen.VariableStruct.Capacity;
qp_form.B.lb = 0;
qp_form.A.H = 0;
qp_form.A.f = 0;
qp_form.B.H = 0;
qp_form.B.f = 0;
end%ends load_ac_dc

function qp_form = load_utility(util,Source)
if any(strcmp(Source, {'NG';'Diesel';}))%if it is a fuel utility
    qp_form = [];
    qp_form.output = [];
    qp_form.states = [];
else
    qp_form.states = {'X'};
    qp_form.X.H = 0;
    qp_form.X.f = 1;
    qp_form.X.lb = util.MinImportThresh;
    qp_form.X.ub = inf;% no sellback allowed (only 1 state)
    if util.MinImportThresh<=0 && (util.SellBackRate>0  || (util.SellBackRate==-1 && util.SellBackPerc~=1)) %add sell back state
        qp_form.states = {'X';'Y'};
        if util.SellBackRate>0 
            qp_form.Y.f = -util.SellBackRate;%constant sell back rate
        else
            qp_form.Y.f = -min(util.SellBackPerc/100,1-1e-6);%ensure less than 1, so no issues with pass through power
        end
        qp_form.Y.H = 0;
        qp_form.Y.lb = 0;
        qp_form.Y.ub = inf;
        qp_form.X.lb = 0;
        if strcmp(Source, 'Electricity')
            qp_form.output.E = [1;-1];
        elseif strcmp(Source, 'Heat')%loads the parameters for a distric heating supply. 
            qp_form.output.H = [1;-1];
        elseif strcmp(Source, 'Cooling')%loads the parameters for a distric cooling supply. 
            qp_form.output.C = [1;-1];
        else 
        end
    else
        if strcmp(Source, 'Electricity')
            qp_form.output.E = 1;
        elseif strcmp(Source, 'Heat')%loads the parameters for a distric heating supply. 
            qp_form.output.H = 1;
        elseif strcmp(Source, 'Cooling')%loads the parameters for a distric cooling supply. 
            qp_form.output.C = 1;
        else 
        end
    end
end
end%Ends function load_utility

function qp_form = load_piecewise(gen,n)
% this function loads the parameters for a combined heat and power
% generator, regular electric generator, chiller, or cooling tower
% n is number of segments
% order is either 1 or 2, when order is 2 it solves for the quadratic 
% coefficients of C = c_0 + a_1*x_1 + a_2*x_2 + ... a_n*X_n + b_1*x_1^2 + b_2*x_2^2 + ... b_n*X_n^2
% subject to b_i >0, and a_i> a_(i-1) + b_(i-1)*(x_i)_max
% if order is 1 it solves for linear coefficients 
% of C = c_0 + a_1*x_1 + a_2*x_2 + ... a_n*X_n
%subject to  a_i> a_(i-1)
switch gen.Type
    case{'CHP Generator';'Electric Generator';'Hydrogen Generator'}
        order = 2;
        if isfield(gen.VariableStruct.Startup,'Electricity')
            lower_bound = gen.VariableStruct.Startup.Electricity(end);
            efficiency = gen.Output.Electricity;
            qp_form.output.E = 1;
        elseif isfield(gen.VariableStruct.Startup,'DirectCurrent')
            lower_bound = gen.VariableStruct.Startup.DirectCurrent(end);
            efficiency = gen.Output.DirectCurrent;
            qp_form.output.DC = 1;
        end
    case 'Electrolyzer'
        order = 2;
        lower_bound = gen.VariableStruct.Startup.Hydrogen(end);
        qp_form.output.Hy = 1;
        efficiency = gen.Output.Hydrogen;
    case 'Heater'
        order = 2;
        lower_bound = gen.VariableStruct.Startup.Heat(end);
        qp_form.output.H = 1;
        efficiency = gen.Output.Heat;
    case 'Chiller'
        order = 1;
        lower_bound = gen.VariableStruct.Startup.Cooling(end);
        qp_form.output.C = 1;
        efficiency = gen.Output.Cooling;
    case 'Cooling Tower'
        order = 1;
        lower_bound = gen.VariableStruct.Startup.heat_reject(end);
        qp_form.output.CW = -1;
        efficiency = gen.Output.heat_reject;
end
upper_bound = gen.Size;
capacity = gen.Output.Capacity;
operation_range = find(capacity>=lower_bound/upper_bound);
[p,i_sort] = sort(capacity(operation_range));
eff = efficiency(operation_range);
y = p./eff(i_sort); %cost of generator in terms of input at outputs P
y(isnan(y)) = 0;

seg_end = lower_bound/upper_bound + (1-lower_bound/upper_bound)*(1:n)/n;
[~,a,~] = fitfcn(p,y,seg_end,n,order,0);% in Fit A: c_0 = 0
[~,b,c_0] = fitfcn(p,y,seg_end,n,order,1);


if strcmp(gen.Type,'CHP Generator')
    [h0,heat_out] = fit_coproduction(gen,a,b,n);
    qp_form.constDemand.H = -h0*upper_bound;
    keep_A = remove_segments(a,n,order,heat_out(:,1));
    keep_B = remove_segments(b,n,order,heat_out(:,2));
    qp_form.output.H = zeros(n,2);
    qp_form.output.H(1:nnz(keep_A(1:n)),1)  = heat_out(keep_A(1:n),1);
    qp_form.output.H(1:nnz(keep_B(1:n)),2)  = heat_out(keep_B(1:n),2);
else
    keep_A = remove_segments(a,n,order,[]);
    keep_B = remove_segments(b,n,order,[]);
end
n_a = nnz(keep_A)/order;
n_b = nnz(keep_B)/order;
if strcmp(gen.Type,'Chiller')
    if strcmp(gen.Source,'Electricity')
        source = 'E';
    else
        source = 'H';
        qp_form.output.E = 0;%absorption chiller has no direct electrical load, but later pump for cooling tower can be added
    end
    qp_form.output.(source)(1:n_a,1) = -a(keep_A);
    qp_form.output.(source)(1:n_b,2) = -b(keep_B);
    qp_form.output.CW(1:n_a,1) = 1+ a(keep_A);
    qp_form.output.CW(1:n_b,2) = 1+ b(keep_B);
    qp_form.constDemand.(source) = c_0*upper_bound;
    qp_form.constDemand.CW = -c_0*upper_bound;%heat added to water loop appears as a negative demand for energy in the water loop
    a = a*0;
    b = b*0;
    c_0 = 0;
elseif strcmp(gen.Type,'Cooling Tower')
    qp_form.output.E(1:nnz(keep_A),1) = -a(keep_A);
    qp_form.output.E(1:nnz(keep_B),2) = -b(keep_B);
    qp_form.constDemand.E = c_0*upper_bound;
    qp_form.constDemand.CW = -c_0*upper_bound;%heat added to water loop appears as a negative demand for energy in the water loop
    a = a*0;
    b = b*0;
    c_0 = 0;
end

%put into form that build matrices will see
letters = {'A';'B';'C';'D';'E';'F';'G';'H';'I';'J';'K';'L';'M';'N';'O';'P';};
qp_form.states = cell(max(n_a,n_b),2);
qp_form.states(1:n_a,1) = letters(1:n_a);
qp_form.states(1:n_b,2) = letters(1:n_b);
qp_form.constCost = c_0*upper_bound;
x_max = seg_end - [0,seg_end(1:n-1)];
j = 0;
for i = 1:1:n
    if keep_A(i)
        j = j+1;
        if order == 2
            qp_form.(letters{j}).H = a(n+i)/upper_bound;%% scale by size of gen (later can normalize by changing QPform.Output so that all states are 0-1)
        else
            qp_form.(letters{j}).H = 0;
        end
        qp_form.(letters{j}).f = a(i);
        qp_form.(letters{j}).lb = 0;
        qp_form.(letters{j}).ub = x_max(i)*upper_bound;
    else
        qp_form.(letters{j}).ub = qp_form.(letters{j}).ub + x_max(i)*upper_bound;
    end
end
j = 0;
for i = 1:1:n
    if keep_B(i)        
        j = j+1;
        if order == 2
            qp_form.(letters{j}).H(2) = b(n+i)/upper_bound;%% scale by size of gen (later can normalize by changing QPform.Output so that all states are 0-1)
        else
            qp_form.(letters{j}).H(2) = 0;
        end
        qp_form.(letters{j}).f(2) = b(i);
        if i ==1
            qp_form.(letters{j}).lb(2) = lower_bound;
            qp_form.(letters{j}).ub(2) = max(lower_bound,x_max(i)*upper_bound);
        else
            qp_form.(letters{j}).lb(2) = 0;
            qp_form.(letters{j}).ub(2) = x_max(i)*upper_bound;
        end
    else
        qp_form.(letters{j}).ub(2) = qp_form.(letters{j}).ub(2) + x_max(i)*upper_bound;
    end
end

% % %% plot
% figure
% plot(P*UB,Y*UB,'b')
% hold on
% X2 = linspace(0,x_max_A(1),10);
% if order == 1
%     Y2 = linspace(0,x_max_A(1),10)*A(1);
% else Y2 = linspace(0,x_max_A(1),10)*A(1) + linspace(0,x_max_A(1),10).^2*A(n_A+1);
% end
% for i = 2:1:n_A
%     X2 = [X2, X2(end)+linspace(0,x_max_A(i),10)];
%     if order == 1
%         Y2 = [Y2,Y2(end)+linspace(0,x_max_A(i),10)*A(i)];
%     else Y2 = [Y2,Y2(end)+linspace(0,x_max_A(i),10)*A(i) + linspace(0,x_max_A(i),10).^2*A(n_A+i)];
%     end
% end
% plot(X2,Y2,'g')
% 
% X3 = linspace(0,x_max_B(1),10);
% if order == 1
%     Y3 = c_0+linspace(0,x_max_B(1),10)*B(1);
% else Y3 = c_0+linspace(0,x_max_B(1),10)*B(1) + linspace(0,x_max_B(1),10).^2*B(n_B+1);
% end
% for i = 2:1:n_B
%     X3 = [X3, X3(end)+linspace(0,x_max_B(i),10)];
%     if order == 1
%         Y3 = [Y3,Y3(end)+linspace(0,x_max_B(i),10)*B(i)];
%     else Y3 = [Y3,Y3(end)+linspace(0,x_max_B(i),10)*B(i) + linspace(0,x_max_B(i),10).^2*B(n_B+i)];
%     end
% end
% plot(X3,Y3,'r')
end%Ends function load_piecewise

function keep = remove_segments(a,n,order,co_prod)
%remove segments that are equal (same slope & quadratic) to previous segment
keep = true(1,order*n);
j = 1;
for i = 2:1:n
    if a(i)-a(i-1)<1e-5 && (order == 1 || a(n+i) < 1e-8) && (isempty(co_prod)==1 || abs(co_prod(i) - co_prod(i-1))<1e-5)
        keep(i) = false;
        if order ==2
            keep(n+i) = false;
        end
    else
        j = j+1;
    end
end
end%ends function remove_segments

function [qp_form,dx_dt] = load_storage(gen,scale)
%this function just directs to either hot or cold thermal storage
%if we can get rid of the CS, HS structures then this function can be
%eliminated, because all storage can be handled the same way.
stor.Size = gen.Size*scale;
stor.SelfDischarge  = gen.VariableStruct.SelfDischarge;% SelfDischarge per hour (fraction of total charge)
if isfield(gen.VariableStruct, 'EnStoreType') %if its thermal/hydrogen storage
    stor.PeakDisch = (gen.VariableStruct.DischRatePerc/100*gen.Size); %Thermal kW out
    stor.PeakCharge = (gen.VariableStruct.FillRatePerc/100*gen.Size); %Thermal kW in
    stor.ChargeEff = gen.VariableStruct.ChargeEff;
    stor.DischEff = gen.VariableStruct.DischargeEff;
    stor.UsableSize  = stor.Size; % usable size 
    if strcmp(gen.VariableStruct.EnStoreType, 'ColdTES')
        qp_form.output.C = 1; 
    elseif strcmp(gen.VariableStruct.EnStoreType, 'HotTES')
        qp_form.output.H = 1;
    elseif strcmp(gen.VariableStruct.EnStoreType, 'Hydrogen')
        qp_form.output.Hy = 1;
    end
    
else %electric battery
    qp_form.output.DC = 1; 
    stor.Voltage = gen.VariableStruct.Voltage;
    DischCurrent = gen.VariableStruct.PeakDisch.*gen.Size/stor.Voltage*1000;
    stor.DischResistScaled = (100/DischCurrent)*gen.VariableStruct.DischResist*(1/1000); %Scale so the loss of power is equivelant to that specified at 100Amps
    DischVoltLoss = DischCurrent*stor.DischResistScaled; %keep in mind when calculating loss as function of discharge current
    ChargeCurrent = gen.VariableStruct.PeakCharge*gen.Size/stor.Voltage*1000;
    stor.ChargeResistScaled = (100/ChargeCurrent)*gen.VariableStruct.ChargeResist*(1/1000); %Scale so the loss of power is equivelant to that specified at 100Amps
    ChargeVoltLoss = ChargeCurrent*stor.ChargeResistScaled;
    stor.ChargeEff = stor.Voltage/(stor.Voltage+ChargeVoltLoss);
    stor.DischEff = (stor.Voltage-DischVoltLoss)/stor.Voltage;
    stor.PeakDisch = (DischCurrent*stor.Voltage*stor.DischEff/1000); %PeakDischPower kW out
    stor.PeakCharge = (ChargeCurrent*stor.Voltage/stor.ChargeEff/1000); % PeakChargePower kW in
    stor.UsableSize  = stor.Size*(gen.VariableStruct.MaxDOD/100); % usable size (state x is 0_> usable size, must calculate this initial condition each time from voltage & SOC
end

qp_form.Stor = stor;
dx_dt = stor.PeakDisch;% storage discharge constraint in kW

qp_form.states = {'X';}; %state of charge, charging power, no buffers
qp_form.X.lb = 0;
qp_form.X.ub = stor.UsableSize;
qp_form.X.H = 0;
qp_form.X.f = 0;
a = (1/stor.ChargeEff - stor.DischEff);
qp_form.Ramp.b = [stor.PeakCharge; stor.PeakDisch];%SOC2-SOC1<PeakCharge; SOC1-SOC2<PeakDisch
if a~=0 %not ideal storage, add charging state
    qp_form.states = {'X'; 'Y'}; %state of charge, charging power, no buffers
    qp_form.link.ineq = [a -1]; % (1/nc - nd)*(SOC(t) -SOC(t-1)) - charging <0 ------ Charging is the 1/inefficiency + 1
    qp_form.link.bineq = 0;
    qp_form.Y.lb = 0;
    qp_form.Y.ub = stor.PeakCharge;%the limit on how much charging power can be delivered is handled by the generators' limits, so put inf here to prevent redundancy
    qp_form.Y.H = 0;%the cost of the charging power is handled by the generators 
    qp_form.Y.f = 0;
end    
if isfield(gen.VariableStruct,'Buffer') && gen.VariableStruct.Buffer ~= 0 %buffer states
    if a==0 %ideal storage, ignore charging state
        qp_form.states = {'X';'U';'L'};%SOC(t+1), charging power, upper buffer, lower buffer
        qp_form.link.ineq = [-1 0 -1; 1 -1 0];%-SOC(t)-lowerbuffer<-0.2  and %SOC-upperbuffer<0.8
        qp_form.link.bineq = [0;0;]; %% note: the magnitude of the buffer is set later in find_buffer
    else
        qp_form.states = {'X';'Y';'U';'L'};%SOC(t+1), charging power, upper buffer, lower buffer
        qp_form.link.ineq = [a, -1, 0, 0];% (1/nc - nd)*(SOC(t) -SOC(t-1)) - charging <0 ------ Charging is the 1/inefficiency + 1
        qp_form.link.ineq = [qp_form.link.ineq; -1 0 0 -1];%-SOC(t)-lowerbuffer<-0.2
        qp_form.link.ineq = [qp_form.link.ineq; 1 0 -1 0];%SOC-upperbuffer<0.8
        qp_form.link.bineq = [0;0;0;];
    end
    qp_form.U.lb = 0;
    qp_form.U.ub = 0;%% note: the magnitude of the buffer is set later in find_buffer
    qp_form.U.H = 0;
    qp_form.U.f = 0;

    qp_form.L.lb = 0;
    qp_form.L.ub = 0;%% note: the magnitude of the buffer is set later in find_buffer
    qp_form.L.H = 0;
    qp_form.L.f = 0;
end
end%Ends function load_storage

function qp_form = load_hydro_storage(gen,scale)
% this function loads the parameters for a hydroelectric plant.
eff = gen.VariableStruct.MaxGenCapacity/(gen.VariableStruct.MaxGenFlow*gen.VariableStruct.MaxHead*84.674);%Power (kW)/ideal power in kW
qp_form.Stor.Size = gen.Size*scale;
qp_form.Stor.SelfDischarge  = 0; %needs to be evaporative losses
qp_form.Stor.DischEff = 1; %100% efficient
qp_form.Stor.UsableSize = qp_form.Stor.Size*gen.VariableStruct.MinHead/gen.VariableStruct.MaxHead;
qp_form.Stor.Power2Flow = 1/(eff*gen.VariableStruct.MaxHead*84.674);%Power (kW) = efficiency(%) * Flow (1000 ft^3/s) * Head (ft) * 84.674 kJ/ (1000ft^3*ft)
qp_form.output.W = 0;
qp_form.output.E = [1;0;];

qp_form.states = {'X';'Y'}; %Power and state of charge
qp_form.X.lb = 0;
qp_form.X.ub = gen.VariableStruct.MaxGenCapacity;
qp_form.X.H = 0;
qp_form.X.f = 0;

qp_form.Y.lb = 0;
qp_form.Y.ub = qp_form.Stor.UsableSize;
qp_form.Y.H = 0;
qp_form.Y.f = 0;

qp_form.Ramp.b = [gen.VariableStruct.RampDown; gen.VariableStruct.RampUp;];%change in power generation

qp_form.link.eq = [qp_form.Stor.Power2Flow, 0];%convert power to flow rate (1000 cfs)
qp_form.link.beq = 0;
    
if gen.VariableStruct.MaxSpillFlow>0
    qp_form.states = {'X';'Y';'S'}; %add spill flow state
    qp_form.S.lb = 0;
    qp_form.S.ub = gen.VariableStruct.MaxSpillFlow;
    qp_form.S.H = 0;
    qp_form.S.f = 0;
    qp_form.link.eq = [qp_form.Stor.Power2Flow, 0, 1];
    qp_form.link.beq = 0;
    qp_form.output.E = [1; 0; 0;];
end
if isfield(gen.VariableStruct,'Buffer') && gen.VariableStruct.Buffer ~= 0 %buffer states
    if gen.VariableStruct.MaxSpillFlow>0
        qp_form.states = {'X';'Y';'S';'U';'L'};
        qp_form.link.eq = [qp_form.Stor.Power2Flow, 0, 1, 0, 0];
        qp_form.link.ineq = [0 -1 0 0  -1; 0 1 0 -1 0;];%-SOC(t)-lowerbuffer<-0.2, %SOC-upperbuffer<0.8
        qp_form.output.E = [1; 0; 0; 0; 0;];
    else
        qp_form.states = {'X';'Y';'U';'L'};
        qp_form.link.eq = [qp_form.Stor.Power2Flow, 0, 0, 0];
        qp_form.link.ineq = [0 -1 0  -1; 0 1 -1 0;];%-SOC(t)-lowerbuffer<-0.2, %SOC-upperbuffer<0.8
        qp_form.output.E = [1; 0; 0; 0;];
    end
    qp_form.link.bineq = [0;0;];

    qp_form.U.lb = 0;
    qp_form.U.ub = 0;
    qp_form.U.H = 0;
    qp_form.U.f = 0;

    qp_form.L.lb = 0;
    qp_form.L.ub = 0;
    qp_form.L.H = 0;
    qp_form.L.f = 0;
end
end%Ends function load_hydro_storage

function [r_square,A,c_0] = fitfcn(output,input,seg_end,n,order,intercept)
x_max = seg_end - [0,seg_end(1:n-1)];
X = zeros(length(output),n);
for i = 1:1:length(output)
    X(i,1) = min(output(i),seg_end(1));
    for j = 2:1:n
        X(i,j) = max(0,min(x_max(j),output(i)-sum(X(i,1:j-1))));
    end
end
QP.H = ones(order*n);
QP.f = ones(order*n,1);
for i = 1:1:order*n
    if i<=n
        a = X(:,i);
    else
        a = X(:,i-n).^2;
    end
    for j = 1:1:order*n
        if j<=n
            b = X(:,j);
        else
            b = X(:,j-n).^2;
        end
        QP.H(i,j) = 2*sum(a.*b);
    end
    QP.f(i) = -2*sum(a.*input);
end

QP.A = zeros(order*n-1,order*n);
QP.b = zeros(order*n-1,1);
for i = 1:1:n-1
    QP.A(i,i) = 1;
    QP.A(i,i+1) = -1;
    if order==2
        QP.A(i,n+i) = 2*x_max(i);
    end
end
if order==2
    for i = n:1:2*n-1
        QP.A(i,i+1) = -1;
    end
    QP.Aeq = [x_max x_max.^2];
else
    QP.Aeq = x_max;
end
QP.beq = input(end);
QP.lb = [];
QP.ub =[];
QP.solver = 'quadprog';
if intercept
    a = zeros(1,order*n);
    for i = 1:1:order*n
        if i<=n
            a(i) = 2*sum(X(:,i));
        else
            a(i) = 2*sum(X(:,i-n).^2);
        end
    end
    QP.H = [2*length(output), a; [a',QP.H];];
    QP.f = [-2*sum(input);QP.f];
    QP.A = [zeros(order*n-1,1),QP.A];
    QP.Aeq = [1 QP.Aeq];
end
A = call_solver(QP);

    
if order == 2
    for i = n+1:2*n
        if A(i) <1e-6
            A(i) = 0;
        end
    end
end
if intercept
    c_0 = A(1);
    A = A(2:end);
else
    c_0 = 0;
end
%find residual of result
SS_res = 0;
SS_mean = 0;
mean_input = mean(input);
for i = 1:1:length(output)
    if order == 2
        fit_i = c_0 + sum(A'.*[X(i,:),X(i,:).^2]);
    else
        fit_i = c_0 + sum(A'.*X(i,:));
    end
    SS_res = SS_res + (input(i) - fit_i)^2;
    SS_mean = SS_mean + (input(i) - mean_input)^2;
end
r_square = 1 - SS_res/SS_mean;
end%Ends fitfcn


function [h0,heat_out] = fit_coproduction(gen,a,b,n)
%sets up two least squares problems to find the coefficients
%of the heat co-production
upper_bound = gen.Size;
capacity = gen.Output.Capacity;
if isfield(gen.VariableStruct.Startup,'Electricity')
    lower_bound = gen.VariableStruct.Startup.Electricity(end);
elseif isfield(gen.VariableStruct.Startup,'DirectCurrent')
    lower_bound = gen.VariableStruct.Startup.DirectCurrent(end);
end
seg_end = lower_bound + (upper_bound-lower_bound)*(1:n)/n;
[cap,I] = sort(capacity);
heat = gen.Output.Heat(I);
if isfield(gen.Output,'Electricity')
    elec = gen.Output.Electricity(I);
elseif isfield(gen.Output,'DirectCurrent')
    elec = gen.Output.DirectCurrent(I);
end
seg_heat = [0,seg_end./interp1(cap,elec,seg_end/upper_bound).*interp1(cap,heat,seg_end/upper_bound)]*0.95;% reduce heat co-production by 5% for fitting, because it is better to slightly underestimate
cumulative_slope = seg_heat(2:end)./seg_end;
seg_heat(2:end) = min(seg_heat(2:end),1.3*mean(cumulative_slope(1:n-1))*seg_end);
local_slope = (seg_heat(2:end)-seg_heat(1:end-1))'./(seg_end-[0,seg_end(1:end-1)])';

max_slope_inc = (a(2:n)-a(1:n-1))./a(1:n-1).*local_slope(1:n-1);%maximum that beta (slope of heat ratio) can increase and remain convex due to convex cost function

seg_5n = [lower_bound*ones(n,1);linspace(lower_bound,upper_bound,5*n+1)']/upper_bound;
seg_heat_5n = seg_5n./interp1(cap,elec,seg_5n).*interp1(cap,heat,seg_5n)*0.95;

%     seg_heat_5n = min(seg_heat_5n,1.3*mean(cumulative_slope(1:n-1))*seg_5n)*0.95;

[out_h1,~] = fitfcn2(seg_5n,seg_heat_5n,seg_end/upper_bound,n,max_slope_inc,[]);

%repeat for non-zero y-intercept
slope_seg1 = (seg_heat(2) - lower_bound./interp1(cap,elec,lower_bound/upper_bound).*interp1(cap,heat,lower_bound/upper_bound))/(seg_end(1)-lower_bound);%slope of 1st segment wih non-zero y-intercept
seg_heat(1) = seg_heat(2) - seg_end(1)*slope_seg1;
local_slope = (seg_heat(2:end)-seg_heat(1:end-1))./(seg_end-[0,seg_end(1:end-1)]);
max_slope_inc = 2*(b(2:n)'-b(1:n-1)')./b(1:n-1)'.*local_slope(1:n-1);%maximum that beta (slope of heat ratio) can increase and remain convex due to convex cost function
[out_h2,h0] = fitfcn2(seg_5n,seg_heat_5n,seg_end/upper_bound,n,max_slope_inc,1.1*local_slope(1));
heat_out = [out_h1,out_h2];
end


function [fit,H_0] = fitfcn2(output,input,seg_end,n,max_slope_inc,local_slope1)
%sets up a least squares problem to find the linear piecewise coefficients
%of the heat co-production
if isempty(local_slope1)
    intercept = false;
else
    intercept = true;
end
x_max = zeros(1,n);
x_max(1) = seg_end(1);
x_max(2:n) = seg_end(2:n) - seg_end(1:n-1);
X = zeros(length(output),n);
for i = 1:1:length(output)
    X(i,1) = min(output(i),seg_end(1));
    for j = 2:1:n
        X(i,j) = max(0,min(x_max(j),output(i)-sum(X(i,1:j-1))));
    end
end
QP.H = ones(n);
QP.f = ones(n,1);
for i = 1:1:n
    a = X(:,i);
    for j = 1:1:n
        b = X(:,j);
        QP.H(i,j) = 2*sum(a.*b);
    end
    QP.f(i) = -2*sum(a.*input);
end

QP.A = zeros(n-1,n);
QP.b = max_slope_inc';
for i = 1:1:n-1 %the sign of this constraint is switched from fitfcn because we want the coefficients of heat recovery to get smaller
    QP.A(i,i) = -1;
    QP.A(i,i+1) = 1;
end
QP.Aeq = [];
QP.beq = [];
QP.lb = [];
QP.ub =[];
QP.solver = 'quadprog';
if intercept
    a = zeros(1,n);
    for i = 1:1:n
        a(i) = 2*sum(X(:,i));
    end
    QP.H = [2*length(X), a; [a',QP.H];];
    QP.f = [-2*sum(input);QP.f];
    QP.A = [zeros(n-1,1),QP.A; 0 1, zeros(1,n-1)];
    QP.b = [QP.b;local_slope1];
end
A = call_solver(QP);
if intercept
    H_0 = A(1);
    fit = A(2:n+1);
else
    H_0 = 0;
    fit = A;
end
end% ends function fitfcn2