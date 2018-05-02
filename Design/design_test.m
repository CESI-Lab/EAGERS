function [npc,costs,mc] = design_test(project,size,index_size,test_data,design_day,years)
%PLANTNPC Calculate a Net Present Cost value for a given set of generators
%specified by iScaleGen, sized as specified by size.
% INPUTS
%   size        Number array. Each generator's new size. Length should
%               equal that of iScaleGen.
%   iScaleGen   Integer array. The indices of those generators that will be
%               scaled. Length should equal that of scale.
%   DesignDay   Boolean. Whether disign day method will be used.
%   Years       Number. Number of years ahead to be considered in NPC
%               calculation.
%
% OUTPUTS
%   npc         Number. The calculated net present cost.

gen = project.Generator;
options = project.optimoptions;
network  = project.Network;
options.method = 'Planning';
options.forecast = 'Perfect';%Perfect forecast pulls directly from TestData
if isfield(project,'Building') && ~isempty(project.Building)
    buildings = project.Building;
    options.forecast = 'Building';
else
    buildings = [];
end
if isfield(project,'cool_tower')
    cool_tower = project.cool_tower;
else
    cool_tower = [];
end
if isfield(project,'Data') 
    data = project.Data;
else
    data = [];
end
equip_costs = project.Costs.Equipment;
discount_rate = project.Costs.DiscountRate;
[gen,equip_costs] = design_resize(gen,equip_costs,size,index_size);


options.Interval = floor(test_data.Timestamp(end)-test_data.Timestamp(1));
freq = 1; %period of repetition (1 = 1 day)
res = options.Resolution/24;
n_o = round(freq/res)+1;
test_data = update_test_data(test_data,data,gen,options);
test_data.RealTimeData = interpolate_data(test_data,options.Resolution*3600,0.00);%create test data at correct frequency
if design_day ==1%If design days option is selected, optimize the 1st day of the month, and assume the rest of the month to be identical
    options.endSOC = 'Initial';%Constrain the final SOC of any storage device to be equal to the initial charge so that days 2-30 of each month do not over depleate storage
    options.Horizon = max(24,options.Horizon);%Make the horizon at least 1 day
    date_1 = test_data.Timestamp(1);%set the starting date
else
    options.endSOC = 'Flexible';%remove constraint on final SOC of storage
end

[gen,buildings,cool_tower,subnet,op_mat_a,op_mat_b,one_step,online] = initialize_optimization(gen,buildings,cool_tower,network,options,test_data.RealTimeData);

if design_day ==1
    forecast_time = date_1+[0;build_time_vector(options)/24];%linspace(Date,DateEnd)';would need to re-do optimization matrices for this time vector
    [num_steps,dispatch,predicted,design,run_data] = pre_allocate_space(gen,buildings,cool_tower,subnet,options,test_data);%create Plant.Design structure with correct space
    test_data.RealTimeData = interpolate_data(test_data,options.Resolution*3600,0.00);%create test data at correct frequency
    STR = 'Optimizing Design Day Dispatch';
    planning_waitbar=waitbar(0,STR,'Visible','on');

    prev_data = get_data(test_data.RealTimeData,linspace((forecast_time(1) - res - freq),forecast_time(1)-res,n_o)',[]);
    now_data = get_data(test_data.RealTimeData,forecast_time(1),[]);
    future_data = get_data(test_data.RealTimeData,forecast_time(2:end),[]);
    [data_t0,gen,~] = update_forecast(gen,buildings,cool_tower,subnet,options,date_1,test_data.HistProf,prev_data,now_data,future_data);
    [gen,cool_tower] = automatic_ic(gen,buildings,cool_tower,subnet,date_1,one_step,options,data_t0);% set the initial conditions
    design.Timestamp(1) = date_1;

    s_i = 1;
    while forecast_time(1)+options.Horizon/24<=test_data.Timestamp(end)%loop to simulate days 1 to n in TestData
        D = datevec(forecast_time(1));
        if s_i == 1 || D(3) == 1  %If it is the first step, or the first of the month run the actual optimization                     
            prev_data = get_data(test_data.RealTimeData,linspace((forecast_time(1) - res - freq),forecast_time(1)-res,n_o)',[]);
            now_data = get_data(test_data.RealTimeData,forecast_time(1),[]);
            future_data = get_data(test_data.RealTimeData,forecast_time(2:end),[]);
            [forecast,gen,buildings] = update_forecast(gen,buildings,cool_tower,subnet,options,forecast_time(2:end),test_data.HistProf,prev_data,now_data,future_data);%% function that creates demand vector with time intervals coresponding to those selected
            Solution = dispatch_loop(gen,buildings,cool_tower,subnet,op_mat_a,op_mat_b,one_step,options,forecast_time,forecast,[]);
        else%otherwise just change the dates and use the previous solution,
            forecast.Timestamp = forecast_time(2:end);
        end   
        [s_i,forecast_time,gen,buildings,cool_tower,design,dispatch,predicted] = dispatch_record(gen,buildings,cool_tower,subnet,options,test_data.RealTimeData,s_i,forecast_time,forecast,Solution,design,dispatch,predicted);%put solution into Plant.Design
        if isfield(subnet,'Hydro')
            n_s = length(forecast_time)-1;
            for n = 1:1:length(subnet.Hydro.nodes)
                test_data.RealTimeData.Hydro.OutFlow(s_i-n_s+1:s_i,n) = design.OutFlow(s_i-n_s+1:s_i,n);
            end
        end
        n_b = length(buildings);
        for i = 1:1:n_b
            buildings(i).Timestamp = 0;
        end
        waitbar(s_i/num_steps,planning_waitbar,strcat('Running Design Day Dispatch'));
    end
    close(planning_waitbar)
else
    STR = 'Optimizing Dispatch Throughout Entire Year';
    planning_waitbar=waitbar(0,STR,'Visible','on');
    [gen,buildings,cool_tower,test_data,design,~,~,~] = run_simulation(test_data.Timestamp(1),0,1,[],test_data.RealTimeData,test_data.HistProf,gen,buildings,cool_tower,network,options,subnet,data);
    close(planning_waitbar)
end
[costs,npc,mc] = design_costs(gen,design,options,years,equip_costs,discount_rate);%update the costs, monthly costs & NPC for system k
end%ends function design_test