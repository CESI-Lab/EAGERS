function [gen,buildings,fluid_loop,design,dispatch,predicted] = run_simulation(date,num_steps,s_i,handles,test_data,hist_prof,gen,buildings,fluid_loop,options,subnet,op_mat_a,op_mat_b,one_step,design,dispatch,predicted)
%num_steps: the number of dispatch optimizations that will occur during the entire simulation
%s_i: Counter for dispatch loop
freq = 1; %period of repetition (1 = 1 day)
res = options.Resolution/24;
n_o = round(freq/res)+1;
time = build_time_vector(options);%% set up vector of time interval
date = date+[0;time/24];
[wy_forecast, gen] = water_year_forecast(gen,buildings,fluid_loop,subnet,options,date,[],test_data);%if october 1st,Run a yearly forecast for hydrology
if num_steps == 0
    s_i = 1;
    [num_steps,dispatch,predicted,design,~] = pre_allocate_space(gen,buildings,fluid_loop,subnet,options,test_data);
    % k = center_menu('Select Option','Manually specify initial conditions','Automatically determine initial conditions');
    k = 2;
    if k ==1
        [gen,buildings] = manual_ic(gen,buildings,fluid_loop,subnet);
    else
        prev_data = get_data(test_data,linspace((date(1) - res - freq),date(1)-res,n_o)',[],[]);
        now_data = get_data(test_data,date(1),[],[]);
        future_data = get_data(test_data,date(2:end),[],[]);
        [data_t0,gen,~] = update_forecast(gen,buildings,fluid_loop,subnet,options,date(1),hist_prof,prev_data,now_data,future_data);
        [gen,fluid_loop] = automatic_ic(gen,buildings,fluid_loop,subnet,date(1),one_step,options,data_t0);% set the initial conditions
    end
    for i = 1:1:length(gen)
        dispatch.GeneratorState(1,i) = gen(i).CurrentState(1);
        if strcmp(gen(i).Type,'Hydro Storage')
            n = gen(i).QPform.Hydro.subnetNode;%dam #
            dispatch.hydroSOC(1,n) = gen(i).CurrentState(2);
        end
    end
    for i = 1:1:length(buildings)
        dispatch.Buildings(1,i) = buildings(i).Tzone;
    end
    dispatch.Timestamp(1) = date(1);
end
timers = zeros(num_steps,3); % To record times set to zeros(1,3), to not record set to [];
solution.Dispatch = [];
sim_waitbar=waitbar(0,'Running Dispatch','Visible','off');

while s_i<num_steps-1
    prev_data = get_data(test_data,linspace((date(1) - res - freq),date(1)-res,n_o)',[],[]);
    now_data = get_data(test_data,date(1),[],[]);
    future_data = get_data(test_data,date(2:end),[],[]);
    [prev_data,now_data] = water_dispatch_history(prev_data,now_data,dispatch,subnet);
    [forecast,gen,buildings] = update_forecast(gen,buildings,fluid_loop,subnet,options,date(2:end),hist_prof,prev_data,now_data,future_data);%% function that creates demand vector with time intervals coresponding to those selected
    [wy_forecast, gen] = water_year_forecast(gen,buildings,fluid_loop,subnet,options,date,wy_forecast,test_data);%if october 1st,Run a yearly forecast for hydrology
    forecast.wy_forecast = wy_forecast;
    [solution,~] = dispatch_loop(gen,buildings,fluid_loop,subnet,op_mat_a,op_mat_b,one_step,options,date,forecast,solution);
    timers(s_i,:) = solution.timers;
    if ~isempty(handles)
        pause(0.001);
        stop = get(handles.Stop,'Value');
        if stop ==1
            return %stop button was pressed
        end
    end    
    
    [c,~,~] = net_cost(gen,solution.Dispatch,date,'Dispatch');
    predicted.Cost(s_i) = sum(c);
    [s_i,date,gen,buildings,fluid_loop,design,dispatch,predicted] = dispatch_record(gen,buildings,fluid_loop,subnet,options,test_data,s_i,date,forecast,solution,design,dispatch,predicted);
    update_gui_status(handles,solution,date-options.Resolution/24,s_i-1,gen,buildings,fluid_loop,subnet,options,dispatch)
    
    waitbar(s_i/num_steps,sim_waitbar,strcat('Running Dispatch'));
end
