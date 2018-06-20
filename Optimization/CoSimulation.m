function CoSimulation(date, handles)
global Plant  TestData 
%% Start EnergyPlus cosimulation
dir = strrep(which('CoSimulation.m'),fullfile('CoSimulation.m'),'');
cd(fullfile(dir,'NRELmethod'))
installMlep
ep = mlepProcess;
ep.arguments = {'B_2ndCoSim_v2', 'USA_VA_Sterling-Washington.Dulles.Intl.AP.724030_TMY3'};
ep.acceptTimeout = 6000;
[status, msg] = ep.start;  
if status ~= 0
    error('Could not start EnergyPlus: %s.', msg);
end
[status, msg] = ep.acceptSocket;
if status ~= 0
    error('Could not connect to EnergyPlus: %s.', msg);
end

%% EAGERS initialize variables
n_g = length(Plant.Generator);%skip initialization
GenStatus = zeros(1,n_g);
lb = zeros(1,n_g);
for i=1:1:n_g
    switch Plant.Generator(i).Type
        case {'Electric Generator';'CHP Generator';}
            lb(i) = Plant.Generator(i).VariableStruct.Startup.Electricity(end);
            GenStatus(i) = 5;
        case 'Chiller'
            lb(i) = Plant.Generator(i).VariableStruct.Startup.Cooling(end);
            GenStatus(i) = 0;
        case 'Heater'
            lb(i) = Plant.Generator(i).VariableStruct.Startup.Heat(end);
            GenStatus(i) = 0;
        case {'Electric Storage';'Thermal Storage';}
            GenStatus(i) = 0.5*Plant.Generator(i).Size*(Plant.Generator(i).VariableStruct.MaxDOD/100); % IC = halfway charged energy storage
    end
    Plant.Generator(i).CurrentState(1) = GenStatus(i);
    Plant.Generator(i).Status = GenStatus(i)>lb(i);
end
n_s = 365*24/Plant.optimoptions.Resolution+1;%# of simulation steps
TestData.RealTimeData.Timestamp = linspace(datenum([2017,1,1]),datenum([2018,1,1]),n_s)';

n_s = round(Plant.optimoptions.Horizon/Plant.optimoptions.Resolution); %number of steps per dispatch
num_steps = Plant.optimoptions.Interval*24+1; %number of simulation steps
n_1hour = round(1/Plant.optimoptions.Resolution);
Plant.Dispatch.Temperature = zeros(num_steps,1);
Plant.Dispatch.Timestamp = zeros(num_steps,1);
Plant.Dispatch.GeneratorState = zeros(num_steps,n_g);
Plant.Predicted.GenDisp = zeros(n_s+1,n_g,num_steps);
Plant.Predicted.Timestamp = zeros(n_s+1,num_steps);
Plant.Predicted.Cost = zeros(num_steps,1);
Plant.Predicted.Demand = [];

if isfield(Plant,'Building') && ~isempty(Plant.Building)
    buildings = Plant.Building;
else
    buildings = [];
end
Plant.Dispatch.GeneratorState(1,:) = GenStatus;
Plant.Dispatch.Timestamp(1) = date;
freq = 1; %period of repetition (1 = 1 day)
res = Plant.optimoptions.Resolution/24;
n_o = round(freq/res)+1;
s_i=1; %counter for # of times dispatch loop has run
Time = build_time_vector(Plant.optimoptions);%% set up vector of time interval
t_EDC = 10;
while s_i<num_steps-1
    %% Read from EnergyPlus
    packet = ep.read;
    if isempty(packet)
        error('Could not read outputs from E+.');
    end
    
    % Parse it to obtain building outputs
    [flag, eptime, outputs] = mlepDecodePacket(packet);
    Plant.Building.Tzone = outputs(1);
    if flag ~= 0, break; end
    if s_i < 24*t_EDC
        SP = [22,22,0];
    else
        %% run optimization
        date = date+[0;Time/24];
        prev_data = get_data(TestData,linspace((date(1) - res - freq),date(1)-res,n_o)',[],[]);
        now_data = get_data(TestData,date(1),[],[]);
        future_data = get_data(TestData,date(2:end),[],[]);
        [forecast,Plant.Generator,buildings] = update_forecast(Plant.Generator,buildings,[],Plant.subNet,Plant.optimoptions,date(2:end),[],prev_data,now_data,future_data);%% function that creates demand vector with time intervals coresponding to those selected
        scale_cost = update_cost(date(2:end),Plant.Generator); %% All feedstock costs were assumed to be 1 when building matrices 
        solution = solver_nrel(Plant.Generator,Plant.optimoptions,GenStatus(1:n_g),Plant.Building.Tzone,forecast,scale_cost);
        solution.Timestamp = date;
        SP = [solution.Buildings.Temperature(n_1hour),solution.Buildings.Temperature(n_1hour),solution.heat_recovery(n_1hour+1)];
    end
    %% Write to inputs of E+
    ep.write(mlepEncodeRealData(2, 0, (s_i-1)*3600, SP));  

    %% Update status (get from E+)
    Plant.Generator(4).CurrentState(1) = outputs(5);%chiller output
    if s_i >= 24*t_EDC
        %% Plot to GUI
        if strcmp(get(handles.uipanelMain1,'Visible'),'on')
            backSteps = min(s_i,Plant.optimoptions.Horizon/1);
            dispatch.Dispatch = Plant.Dispatch.GeneratorState(s_i-backSteps+1:s_i,:);
            dispatch.Timestamp = Plant.Dispatch.Timestamp(s_i-backSteps+1:s_i,:);
            dispatch.Buildings = Plant.Dispatch.Temperature(s_i-backSteps+1:s_i,:);
            plot_project(Plant.Generator,Plant.Building,[],3,dispatch,solution,[],handles,'Result')
            pause(0.01)
        elseif strcmp(get(handles.uipanelMain2,'Visible'),'on')
            Plant.Market.MarginCost = marginal_cost(Plant.Generator,dispatch,date);
            plotMarginalCapacityCost(handles)
        end

        %% Record data
        Plant.Predicted.GenDisp(:,:,s_i) = solution.Dispatch;
        Plant.Predicted.Timestamp(:,s_i) = date;
        
        Plant.Dispatch.GeneratorState(s_i+1,:) = sum(solution.Dispatch(2:n_1hour,:),1)/n_1hour;
        GenStatus = solution.Dispatch(n_1hour+1,:);
    else
        Plant.Dispatch.GeneratorState(s_i+1,:) = GenStatus;
    end
    s_i = s_i+1;    
    date = round(1e5*(date(1)+1/24))/1e5;%% count forward 1 hour, rounded to nearest second
    Plant.Dispatch.Timestamp(s_i) = date(1);
    Plant.Dispatch.Temperature(s_i,:) = outputs(1);%solution.Buildings.Temperature(n_1hour,:);
end
% Stop EnergyPlus
ep.stop;
disp(['Stopped with flag ' num2str(flag)]);
% ==========FLAGS==============
% Flag	Description
% +1	Simulation reached end time.
% 0	    Normal operation.
% -1	Simulation terminated due to an unspecified error.
% -10	Simulation terminated due to an error during the initialization.
% -20	Simulation terminated due to an error during the time integration.
end%Ends function CoSimulation