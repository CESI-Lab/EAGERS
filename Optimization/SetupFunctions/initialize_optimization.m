function [gen,buildings,fluid_loop,subnet,op_mat_a,op_mat_b,one_step,online,num_steps,dispatch,predicted,design,flag] = initialize_optimization(gen,buildings,fluid_loop,network,options,test_data,date)
%% Load generators, & build QP matrices
Time = build_time_vector(options);%% set up dt vector of time interval length
dt = Time - [0; Time(1:end-1)];
[gen,~] = check_ac_dc(gen,buildings,test_data);%just in case planning tool didn't check this
[gen, buildings] = update_qpform_all(gen,buildings,network,options.scaletime);% updates the QPform field in all generators and buildings
[subnet,gen,buildings] = load_network(network,gen,buildings);
gen = max_utility_sellback(gen,subnet,test_data);
gen = find_buffer(gen,subnet,options.Horizon);%need to wait until everything else has been loaded to load the buffer, because it is reliant on how quickly everything else can charge the system

op_mat_a = load_matrices(gen,buildings,fluid_loop,subnet,options,'A',dt); %build quadratic programming matrices for FitA
op_mat_b = load_matrices(gen,buildings,fluid_loop,subnet,options,'B',dt);%build quadratic programming matrices for FitB
one_step = load_matrices(gen,buildings,fluid_loop,subnet,options,'B',[]);%build quadratic programming matrices for single time step
online = [];
if strcmp(options.method,'Control')
    A.Horizon = options.Resolution;%the horizon is the resolution
    A.Resolution = options.Topt/3600;%the resolution is the frequency of Topt
    A.tspacing = 'constant';
    OnlineTime = build_time_vector(A);%% set up dt vector of time interval length
    dt2 = OnlineTime - [0, OnlineTime(1:end-1)];
    for t = 1:1:length(OnlineTime)
        online(t) = load_matrices(gen,buildings,fluid_loop,subnet,options,'B',dt2(t:end)); %build the matrix for the onlineOptimLoop using FitB
    end
end

[num_steps,dispatch,predicted,design,~] = pre_allocate_space(gen,buildings,fluid_loop,subnet,options,test_data);
% k = center_menu('Select Option','Manually specify initial conditions','Automatically determine initial conditions');
k = 2;
if k ==1
    [gen,buildings] = manual_ic(gen,buildings,fluid_loop,subnet);
    flag = 0;
else
    [data_t0,gen,~] = update_forecast(gen,buildings,fluid_loop,subnet,options,date,test_data,[]);
    [gen,fluid_loop,flag] = automatic_ic(gen,buildings,fluid_loop,subnet,date,one_step,options,data_t0);% set the initial conditions
end

for i = 1:1:length(buildings)
    dispatch.Buildings(1,i) = buildings(i).Tzone;
end
dispatch.Timestamp(1) = date;
end%ends function initialize_optimization
