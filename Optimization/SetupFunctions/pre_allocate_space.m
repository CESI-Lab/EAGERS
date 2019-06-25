function [num_steps,dispatch,predicted,design,run_data] = pre_allocate_space(gen,building,fluid_loop,subnet,options,test_data)
n_g = length(gen);
n_b = length(building);
switch options.method
    case 'Dispatch'
        mode = 'Dispatch';
    case 'Planning'
        mode = 'Design';
    case 'Control'
        mode = 'Control';
end

network_names = fieldnames(subnet);
n_l = 0;
for net = 1:1:length(network_names)
    n_l = n_l + length(subnet.(network_names{net}).lineNames);
end
n_d = 0;
for i = 1:1:n_g
    if strcmp(gen(i).Type,'Hydro Storage')
        n_d = n_d + 1;
    end
end
n_fl = length(fluid_loop);
n_s = round(options.Horizon/options.Resolution); %number of steps per dispatch
num_steps = options.Interval*24/options.Resolution+1; %number of simulation steps
dispatch = [];
predicted = [];
design = [];
run_data = [];
f = fieldnames(test_data);
f = f(~strcmp('Timestamp',f));
if strcmp(mode,'Design')
    design.Timestamp = zeros(num_steps,1);   
    design.GeneratorState = zeros(num_steps,n_g);
    design.LineFlows = zeros(num_steps,n_l);
    design.Buildings = zeros(num_steps,n_b);
    design.fluid_loop = zeros(num_steps,n_fl);
    design.Hydro.SOC = zeros(num_steps,n_d);
    design.Hydro.spillFlow = zeros(num_steps,n_d);
    design.Hydro.OutFlow = zeros(num_steps,n_d);
    design.LBRelax = ones(num_steps,1);%degree to which you relax the lower bound when doing cQP
    for j = 1:1:length(f)
        if isstruct(test_data.(f{j}))
            S = fieldnames(test_data.(f{j}));
            for i = 1:1:length(S)
                if isnumeric(test_data.(f{j}).(S{i}))
                    l = length(test_data.(f{j}).(S{i})(1,:));
                    design.(f{j}).(S{i}) = zeros(num_steps,l);
                end
            end
        elseif ~isempty(test_data.(f{j})) && isnumeric(test_data.(f{j}))
            l = length(test_data.(f{j})(1,:));
            design.(f{j}) = zeros(num_steps,l);
        end
    end
elseif strcmp(mode,'Dispatch')
    dispatch.Timestamp = zeros(num_steps,1);
    dispatch.GeneratorState = zeros(num_steps,n_g);
    dispatch.LineFlows = zeros(num_steps,n_l);
    dispatch.Buildings = zeros(num_steps,n_b);
    dispatch.fluid_loop = zeros(num_steps,n_fl);
    dispatch.Hydro.SOC = zeros(num_steps,n_d);
    dispatch.Hydro.spillFlow = zeros(num_steps,n_d);
    dispatch.Hydro.OutFlow = zeros(num_steps,n_d);
    predicted.GenDisp = zeros(n_s,n_g,num_steps);
    predicted.LineFlows = zeros(n_s,n_l,num_steps);
    predicted.Buildings = zeros(n_s,n_b,num_steps);
    predicted.fluid_loop = zeros(n_s,n_fl,num_steps);
    predicted.hydroSOC = zeros(n_s,n_d,num_steps);
    predicted.Timestamp = zeros(n_s,num_steps);
    predicted.Cost = zeros(num_steps,1);
    predicted.LBRelax = ones(num_steps,1);%degree to which you relax the lower bound when doing cQP
    for j = 1:1:length(f)
        if isstruct(test_data.(f{j}))
            S = fieldnames(test_data.(f{j}));
            for i = 1:1:length(S)
                if isnumeric(test_data.(f{j}).(S{i}))
                    l = length(test_data.(f{j}).(S{i})(1,:));
                    predicted.(f{j}).(S{i}) = zeros(num_steps,n_s,l);
                    dispatch.(f{j}).(S{i}) = zeros(num_steps,l);
                end
            end
        elseif ~isempty(test_data.(f{j})) && isnumeric(test_data.(f{j}))
            l = length(test_data.(f{j})(1,:));
            predicted.(f{j}) = zeros(num_steps,n_s,l);
            dispatch.(f{j}) = zeros(num_steps,l);
        end
    end
    run_data = dispatch;
end
end%Ends pre_allocate_space
