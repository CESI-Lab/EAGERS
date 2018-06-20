function [test_sys,timestamp,costs,npc,mc] = run_planning(test_sys,test_data,years,design_day,base)
for i_ts = 1:1:length(test_sys)%Run through list of projects
    gen = test_sys(i_ts).Generator;
    options = test_sys(i_ts).optimoptions;
    if ~isfield(test_sys(i_ts),'Design') || ~isfield(test_sys(i_ts).Design,'Timestamp') || any(test_sys(i_ts).Design.Timestamp==0)%if the project has already been run, don't re-simulate (need to empty Design when something in the GUI changes what the solution should be
        network  = test_sys(i_ts).Network;
        options.method = 'Planning';
        options.forecast = 'Perfect';%Perfect forecast pulls directly from TestData
        if isfield(test_sys(i_ts),'Building') && ~isempty(test_sys(i_ts).Building)
            buildings = test_sys(i_ts).Building;
            options.forecast = 'Building';
        else
            buildings = [];
        end
        if isfield(test_sys(i_ts),'cool_tower')
            cool_tower = test_sys(i_ts).cool_tower;
        else
            cool_tower = [];
        end
        if isfield(test_sys(i_ts),'Data') 
            data = test_sys(i_ts).Data;
        else
            data = [];
        end
        options.Interval = floor(test_data.Timestamp(end)-test_data.Timestamp(1));
        freq = 1; %period of repetition (1 = 1 day)
        res = options.Resolution/24;
        n_o = round(freq/res)+1;
        test_data = update_test_data(test_data,data,gen,options);
        if design_day ==1%If design days option is selected, optimize the 1st day of the month, and assume the rest of the month to be identical
            options.endSOC = 'Initial';%Constrain the final SOC of any storage device to be equal to the initial charge so that days 2-30 of each month do not over depleate storage
            options.Horizon = max(24,options.Horizon);%Make the horizon at least 1 day
            date_1 = test_data.Timestamp(1);%set the starting date
        else
            options.endSOC = 'Flexible';%remove constraint on final SOC of storage
        end

        [gen,buildings,cool_tower,subnet,op_mat_a,op_mat_b,one_step,online] = initialize_optimization(gen,buildings,cool_tower,network,options,test_data);

        if design_day ==1
            date = date_1+[0;build_time_vector(options)/24];%linspace(Date,DateEnd)';would need to re-do optimization matrices for this time vector
            [wy_forecast, gen] = water_year_forecast(gen,buildings,cool_tower,subnet,options,date,[],test_data);%if october 1st,Run a yearly forecast for hydrology
            [num_steps,dispatch,predicted,design,run_data] = pre_allocate_space(gen,buildings,cool_tower,subnet,options,test_data);%create Plant.Design structure with correct space
            STR = 'Optimizing Design Day Dispatch';
            planning_waitbar=waitbar(0,STR,'Visible','on');
            
            prev_data = get_data(test_data,linspace((date(1) - res - freq),date(1)-res,n_o)',[],[]);
            now_data = get_data(test_data,date(1),[],[]);
            future_data = get_data(test_data,date(2:end),[],[]);
            [data_t0,gen,~] = update_forecast(gen,buildings,cool_tower,subnet,options,date_1,test_data.HistProf,prev_data,now_data,future_data);
            [gen,cool_tower] = automatic_ic(gen,buildings,cool_tower,subnet,date_1,one_step,options,data_t0);% set the initial conditions
            design.Timestamp(1) = date_1;
            
            s_i = 1;
            while date(1)+options.Horizon/24<=test_data.Timestamp(end)%loop to simulate days 1 to n in TestData
                D = datevec(date(1));
                if s_i == 1 || D(3) == 1  %If it is the first step, or the first of the month run the actual optimization                     
                    prev_data = get_data(test_data,linspace((date(1) - res - freq),date(1)-res,n_o)',[],[]);
                    now_data = get_data(test_data,date(1),[],[]);
                    future_data = get_data(test_data,date(2:end),[],[]);
                    [forecast,gen,buildings] = update_forecast(gen,buildings,cool_tower,subnet,options,date(2:end),test_data.HistProf,prev_data,now_data,future_data);%% function that creates demand vector with time intervals coresponding to those selected
                    [wy_forecast, gen] = water_year_forecast(gen,buildings,cool_tower,subnet,options,date,wy_forecast,test_data);%if october 1st,Run a yearly forecast for hydrology
                    forecast.wy_forecast = wy_forecast;
                    [Solution,flag] = dispatch_loop(gen,buildings,cool_tower,subnet,op_mat_a,op_mat_b,one_step,options,date,forecast,[]);
                else%otherwise just change the dates and use the previous solution,
                    forecast.Timestamp = date(2:end);
                end   
                [s_i,date,gen,buildings,cool_tower,design,dispatch,predicted] = dispatch_record(gen,buildings,cool_tower,subnet,options,test_data,s_i,date,forecast,Solution,design,dispatch,predicted);%put solution into Plant.Design
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
            close(planning_waitbar)
        else
            if ~isfield(test_sys(i_ts),'Design') || isempty(test_sys(i_ts).Design) || any(test_sys(i_ts).Design.Timestamp==0) %at least some points have not been run
                STR = 'Optimizing Dispatch Throughout Entire Year';
                planning_waitbar=waitbar(0,STR,'Visible','on');
                [gen,buildings,cool_tower,test_data,design,~,~,~] = run_simulation(test_data.Timestamp(1),0,1,[],test_data,test_data.HistProf,gen,buildings,cool_tower,network,options,subnet,data);
                close(planning_waitbar)
            end
        end
        options.method = test_sys(i_ts).optimoptions.method;%change back method so it is correct when switching to control tool
        options.forecast = test_sys(i_ts).optimoptions.forecast;%change back method so it is correct when switching to control tool
        test_sys(i_ts).Design = design;
        test_sys(i_ts).Generator = gen;
        test_sys(i_ts).optimoptions = options;
        test_sys(i_ts).Building = buildings;
        test_sys(i_ts).cool_tower = cool_tower;
        timestamp(:,i_ts) = design.Timestamp;
    else
        design = test_sys(i_ts).Design;
        timestamp(:,i_ts) = design.Timestamp;
        flag = 0;
    end
    [costs(:,:,i_ts),npc(:,i_ts),mc(:,i_ts)] = design_costs(gen,design,options,years,test_sys(i_ts).Costs.Equipment,test_sys(i_ts).Costs.DiscountRate,flag); % Update the costs, monthly costs & NPC for system k.
end

for i_ts = 1:1:length(test_sys)%calculate caparative economic terms
    if i_ts == base
        test_sys(i_ts).Costs.Financial.irr = 0;
        test_sys(i_ts).Costs.Financial.payback = 0;
    else
        [test_sys(i_ts).Costs.Financial.irr,test_sys(i_ts).Costs.Financial.payback] = financial_metric(test_sys(i_ts).Costs.ProjectedMonthlyCosts,test_sys(base).Costs.ProjectedMonthlyCosts,(1+test_sys(i_ts).Costs.DiscountRate/100));
    end
end
end%ends function run_planning
