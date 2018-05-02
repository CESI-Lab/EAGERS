function [gen,buildings,cool_tower,subnet,op_mat_a,op_mat_b,one_step,online] = initialize_optimization(gen,buildings,cool_tower,network,options,test_data)
%% Load generators, & build QP matrices
Time = build_time_vector(options);%% set up dt vector of time interval length
dt = Time - [0; Time(1:end-1)];
if isfield(test_data,'Demand')
    demand_types = fieldnames(test_data.Demand);
else
    demand_types = {};
end
[gen,~] = check_ac_dc(gen,buildings,demand_types);%just in case planning tool didn't check this
[gen, buildings] = update_qpform_all(gen,buildings,network,options.scaletime);% updates the QPform field in all generators and buildings
[subnet, gen] = load_network(network,gen);
if isfield(network,'Hydro')
    [gen,subnet] = locate_dams(gen,subnet,dt);
end
if ~isempty(buildings)
    buildings = locate_buildings(buildings,subnet,gen);
end
if isfield(test_data,'Demand')
    d_fields = fieldnames(test_data.Demand);
    for i = 1:1: length(d_fields)
        d_max.(d_fields{i}) = max(test_data.Demand.(d_fields{i}));
    end
else
    d_fields = {};
    d_max = [];
end
gen = max_utility_sellback(gen,subnet,d_fields,d_max);
gen = find_buffer(gen,subnet,options.Horizon);%need to wait until everything else has been loaded to load the buffer, because it is reliant on how quickly everything else can charge the system

op_mat_a = load_matrices(gen,buildings,cool_tower,subnet,options,'A',dt); %build quadratic programming matrices for FitA
op_mat_b = load_matrices(gen,buildings,cool_tower,subnet,options,'B',dt);%build quadratic programming matrices for FitB
one_step = load_matrices(gen,buildings,cool_tower,subnet,options,'B',[]);%build quadratic programming matrices for single time step
online = [];
if strcmp(options.method,'Control')
    A.Horizon = options.Resolution;%the horizon is the resolution
    A.Resolution = options.Topt/3600;%the resolution is the frequency of Topt
    A.tspacing = 'constant';
    OnlineTime = build_time_vector(A);%% set up dt vector of time interval length
    dt2 = OnlineTime - [0, OnlineTime(1:end-1)];
    for t = 1:1:length(OnlineTime)
        online(t) = load_matrices(gen,buildings,cool_tower,subnet,options,'B',dt2(t:end)); %build the matrix for the onlineOptimLoop using FitB
    end
end
end%ends function initialize_optimization
