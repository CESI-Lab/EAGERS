function [gen,buildings,dispatch,predicted] = run_simulation_NREL(gen,buildings,test_data,options,date,handles,ep)
%% initialize variables
n_g = length(gen);%skip initialization
lb = zeros(1,n_g);
for i=1:1:n_g
    switch gen(i).Type
        case {'Electric Generator';'CHP Generator';}
            lb(i) = gen(i).VariableStruct.Startup.Electricity(end);
            gen(i).CurrentState(1) = 5;
        case 'Chiller'
            lb(i) = gen(i).VariableStruct.Startup.Cooling(end);
            gen(i).CurrentState(1) = 0;
        case 'Heater'
            lb(i) = gen(i).VariableStruct.Startup.Heat(end);
            gen(i).CurrentState(1) = 0;
        case {'Electric Storage';'Thermal Storage';}
            gen(i).CurrentState(1) = 0.5*gen(i).Size*(gen(i).VariableStruct.MaxDOD/100); % IC = halfway charged energy storage
        otherwise
            gen(i).CurrentState(1) = 0;
    end
    gen(i).Status = gen(i).CurrentState(1)>lb(i);
end

n_s = round(options.Horizon/options.Resolution); %number of steps per dispatch
num_steps = options.Interval*24+1; %number of simulation steps
n_1hour = round(1/options.Resolution);
dispatch.Temperature = zeros(num_steps,1);
dispatch.Timestamp = zeros(num_steps,1);
dispatch.GeneratorState = zeros(num_steps,n_g);
predicted.GenDisp = zeros(n_s+1,n_g,num_steps);
predicted.Timestamp = zeros(n_s+1,num_steps);
predicted.Cost = zeros(num_steps,1);
predicted.Demand = [];

for i = 1:1:n_g
    dispatch.GeneratorState(1,i) = gen(i).CurrentState(1);
end
dispatch.Timestamp(1) = date;
s_i=1; %counter for # of times dispatch loop has run
t_vec = build_time_vector(options);%% set up vector of time interval
t_EDC = 10;
while s_i<num_steps-1
    if options.EnergyPlus
        %% Read from EnergyPlus
        packet = ep.read;
        if isempty(packet)
            error('Could not read outputs from E+.');
        end
        [flag, eptime, outputs] = mlepDecodePacket(packet);
        gen(4).CurrentState(1) = outputs(5);%chiller output
        buildings.Tzone = outputs(1);
        if flag ~= 0, break; end
    end
    if s_i >= 24*t_EDC
        %% run optimization
        date = date+[0;t_vec/24];
        [forecast,gen,buildings] = update_forecast(gen,buildings,[],[],options,date(2:end),test_data,[]);%% function that creates demand vector with time intervals coresponding to those selected
        scale_cost = update_cost(date(2:end),gen); %% All feedstock costs were assumed to be 1 when building matrices 
        solution = solver_nrel(gen,options,buildings.Tzone,forecast,scale_cost);
        solution.Timestamp = date;
    end
    if options.EnergyPlus
        %% Write to inputs of E+
        if s_i < 24*t_EDC
            SP = [22,22,0];
        else
            SP = [solution.Buildings.Temperature(n_1hour),solution.Buildings.Temperature(n_1hour),solution.heat_recovery(n_1hour+1)];
        end
        ep.write(mlepEncodeRealData(2, 0, (s_i-1)*3600, SP));  
    end

    if s_i >= 24*t_EDC
        %% Plot to GUI
        if strcmp(get(handles.uipanelMain1,'Visible'),'on')
            backSteps = min(s_i,options.Horizon/1);
            history.Dispatch = dispatch.GeneratorState(s_i-backSteps+1:s_i,:);
            history.Timestamp = dispatch.Timestamp(s_i-backSteps+1:s_i,:);
            history.Buildings = dispatch.Temperature(s_i-backSteps+1:s_i,:);
            plot_project(gen,buildings,[],3,history,solution,[],handles,'Result')
            pause(0.01)
        end

        %% Record data
        predicted.GenDisp(:,:,s_i) = solution.Dispatch;
        predicted.Timestamp(:,s_i) = date;
        
        dispatch.GeneratorState(s_i+1,:) = sum(solution.Dispatch(2:n_1hour,:),1)/n_1hour;
        for i = 1:1:n_g
        	gen(i).CurrentState(1) = solution.Dispatch(n_1hour+1,i);
        end
    else
        dispatch.GeneratorState(s_i+1,:) = dispatch.GeneratorState(s_i,:);
    end
    s_i = s_i+1;    
    date = round(1e5*(date(1)+1/24))/1e5;%% count forward 1 hour, rounded to nearest second
    dispatch.Timestamp(s_i) = date(1);
    dispatch.Temperature(s_i,:) = outputs(1);%solution.Buildings.Temperature(n_1hour,:);
    if ~isempty(handles)
        pause(0.001);
        stop = get(handles.Stop,'Value');
        if stop ==1
            ep.stop;
            disp(['Stopped with flag ' num2str(flag)]);
            return %stop button was pressed
        end
    end  
end
if options.EnergyPlus
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
end
end%Ends function run_simulation_NREL