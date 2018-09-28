function [npc, costs, mc] = design_test(project, size, multiple, index_gen, test_data, design_day, years)
% Calculate a Net Present Cost value for a given set of generators
% INPUTS
%   size        Number array. Each generator's new size. Length should
%               equal that of index_size.
%   multiple    Number array. how many of each generator system to include.
%   index_gen  Integer array. The indices of those generators that will be
%               scaled. Length should equal that of scale.
%   design_day   Boolean. Whether disign day method will be used.
%   years       Number. Number of years ahead to be considered in NPC
%               calculation.
%
% OUTPUTS
%   npc         Number. The calculated net present cost.
%   costs       array of monthly costs for 1 year broken out by category  (financing, maintinence, etc.)
%   mc          vector of total monthly expenditures for the years of study

%%re-size equipment
[gen,equip_costs,network] = design_resize(project.Generator,project.Costs.Equipment,project.Network,size,multiple,index_gen);

%break the project structure into the variables needed for a single case
options = project.optimoptions;
options.method = 'Planning';%makes the simulations step forward an entire day instead of 1 hour after optimizing the day.
options.forecast = 'Perfect';%Perfect forecast pulls directly from test_data
options.Interval = floor(test_data.Timestamp(end)-test_data.Timestamp(1));
if isfield(project,'Building') && ~isempty(project.Building)
    buildings = project.Building;
else
    buildings = [];
end
if isfield(project,'fluid_loop')
    fluid_loop = project.fluid_loop;
else
    fluid_loop = [];
end
data = [];
if isfield(project,'Data')
    data = project.Data;
end
weather = [];
if isfield(project,'Data')
    weather = project.Data;
end

test_data = update_test_data(test_data,data,weather,options);
if design_day ==1%If design days option is selected, optimize the 1st day of the month, and assume the rest of the month to be identical
    options.endSOC = 'Initial';%Constrain the final SOC of any storage device to be equal to the initial charge so that days 2-30 of each month do not over charge or depleate storage
    options.Horizon = max(24,options.Horizon);%Make the horizon at least 1 day
    date_1 = test_data.Timestamp(1);%set the starting date
else
    options.endSOC = 'Flexible';%remove constraint on final SOC of storage
end
s_i = 1;
[gen,buildings,fluid_loop,subnet,op_mat_a,op_mat_b,one_step,online,num_steps,dispatch,predicted,design,flag] = initialize_optimization(gen,buildings,fluid_loop,network,options,test_data,date_1);
if flag>0
   return
end
if design_day ==1
    date = date_1+[0;build_time_vector(options)/24];%linspace(Date,DateEnd)';would need to re-do optimization matrices for this time vector
    [wy_forecast, gen] = water_year_forecast(gen,buildings,fluid_loop,subnet,options,date,[],test_data);%if october 1st,Run a yearly forecast for hydrology
    STR = 'Optimizing Design Day Dispatch';
    planning_waitbar=waitbar(0,STR,'Visible','on');
    design.Timestamp(1) = date_1;

    while date(1)+options.Horizon/24<=test_data.Timestamp(end)%loop to simulate days 1 to n in TestData
        D = datevec(date(1));
        if s_i == 1 || D(3) == 1
            % It is the first step, or the first of the month; run the
            % actual optimization.
            [forecast,gen,buildings] = update_forecast(gen,buildings,fluid_loop,subnet,options,date(2:end),test_data,dispatch);%% function that creates demand vector with time intervals coresponding to those selected
            [wy_forecast, gen] = water_year_forecast(gen,buildings,fluid_loop,subnet,options,date,wy_forecast,test_data);%if october 1st,Run a yearly forecast for hydrology
            forecast.wy_forecast = wy_forecast;
            [Solution, flag] = dispatch_loop(gen,buildings,fluid_loop,subnet,op_mat_a,op_mat_b,one_step,options,date,forecast,[]);
        else
            % Just change the dates and use the previous solution.
            forecast.Timestamp = date(2:end);
        end   
        if flag>0
            break
        else
            [s_i,date,gen,buildings,fluid_loop,design,dispatch,predicted] = dispatch_record(gen,buildings,fluid_loop,subnet,options,test_data,s_i,date,forecast,Solution,design,dispatch,predicted);%put solution into Plant.Design
            if isfield(subnet,'Hydro')
                n_s = length(date)-1;
                for n = 1:1:length(subnet.Hydro.nodes)
                    test_data.Hydro.OutFlow(s_i-n_s+1:s_i,n) = design.OutFlow(s_i-n_s+1:s_i,n);
                end
            end
            n_b = length(buildings);
            for i = 1:1:n_b
                buildings(i).Timestamp = 0;
            end
            waitbar(s_i/num_steps,planning_waitbar,strcat('Running Design Day Dispatch'));
        end
    end
    close(planning_waitbar)
else
    STR = 'Optimizing Dispatch Throughout Entire Year';
    planning_waitbar=waitbar(0,STR,'Visible','on');
    [gen,buildings,fluid_loop,test_data,design,~,~,~] = run_simulation(test_data.Timestamp(1),num_steps,s_i,[],test_data,test_data.HistProf,gen,buildings,fluid_loop,network,options,subnet,data,[]);
    close(planning_waitbar)
end
[costs,npc,mc] = design_costs(gen,design,options,years,equip_costs,project.Costs.DiscountRate,flag); % Update the costs, monthly costs & NPC for system k.
end % Ends function design_test
