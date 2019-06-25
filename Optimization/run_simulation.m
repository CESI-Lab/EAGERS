function [gen,buildings,fluid_loop,design,dispatch,predicted] = run_simulation(date,num_steps,s_i,handles,test_data,gen,buildings,fluid_loop,options,subnet,op_mat_a,op_mat_b,one_step,design,dispatch,predicted,ep)
%num_steps: the number of dispatch optimizations that will occur during the entire simulation
%s_i: Counter for dispatch loop
time = build_time_vector(options);%% set up vector of time interval
date = date+[0;time/24];
[wy_forecast, gen] = water_year_forecast(gen,buildings,fluid_loop,subnet,options,date,[],test_data);%if october 1st,Run a yearly forecast for hydrology
for i = 1:1:length(gen)
    dispatch.GeneratorState(1,i) = gen(i).CurrentState(1);
    if strcmp(gen(i).Type,'Hydro Storage')
        n = gen(i).QPform.Hydro.subnetNode;%dam #
        dispatch.Hydro.SOC(1,n) = gen(i).CurrentState(2);
    end
end
timers = zeros(num_steps,3); % To record times set to zeros(1,3), to not record set to [];
solution.Dispatch = [];
sim_waitbar=waitbar(0,'Running Dispatch','Visible','off');

while s_i<num_steps-1
    if options.EnergyPlus
        %% Read from EnergyPlus
        packet = ep.read;
        if isempty(packet)
            error('Could not read outputs from E+.');
        end
        [flag, eptime, outputs] = mlepDecodePacket(packet);
        %% need to automatically connect Eplus outputs to status
        buildings.Tzone = outputs(1);
        if flag ~= 0, break; end
    end
    [forecast,gen,buildings] = update_forecast(gen,buildings,fluid_loop,subnet,options,date(2:end),test_data,dispatch);%% function that creates demand vector with time intervals coresponding to those selected
    [wy_forecast, gen] = water_year_forecast(gen,buildings,fluid_loop,subnet,options,date,wy_forecast,test_data);%if october 1st,Run a yearly forecast for hydrology
    forecast.wy_forecast = wy_forecast;
    [solution,~] = dispatch_loop(gen,buildings,fluid_loop,subnet,op_mat_a,op_mat_b,one_step,options,date,forecast,solution);
    timers(s_i,:) = solution.timers;
    if ~isempty(handles)
        pause(0.001);
        stop = get(handles.Stop,'Value');
        if stop ==1
            if options.EnergyPlus
                % Stop EnergyPlus
                ep.stop;
            end
            return %stop button was pressed
        end
    end    
    if options.EnergyPlus
        %% Write to inputs of E+
        if s_i < 24*t_EDC
            SP = [22,22,0];
        else
            n = round(1/options.Resolution);
            heat_recovery = chp_heat(gen(3),solution.Dispatch(1+n,3));
            SP = [solution.Buildings.Temperature(n),solution.Buildings.Temperature(n),heat_recovery];
        end
        ep.write(mlepEncodeRealData(2, 0, (s_i-1)*3600, SP));  
    end
    [c,~,~] = net_cost(gen,solution.Dispatch,date,'Dispatch');
    predicted.Cost(s_i) = sum(c);
    [s_i,date,gen,buildings,fluid_loop,design,dispatch,predicted] = dispatch_record(gen,buildings,fluid_loop,subnet,options,test_data,s_i,date,forecast,solution,design,dispatch,predicted);
    update_gui_status(handles,solution,date-options.Resolution/24,s_i-1,gen,buildings,fluid_loop,subnet,options,dispatch)
    
    waitbar(s_i/num_steps,sim_waitbar,strcat('Running Dispatch'));
end
if options.EnergyPlus% Stop EnergyPlus
    ep.stop;
end