function [si, date, gen, building, fluid_loop, design, dispatch, predicted] ...
    = dispatch_record(gen, building, fluid_loop, subnet, options, test_data, si, date, forecast, solution, design, dispatch, predicted)
% Records the last dispatch and steps forward by a single step or a full day
n_s = length(date)-1;
n_g = length(gen);
n_b = length(building);
n_fl = length(fluid_loop);
dt = (24*3600) * (date(2:end) - date(1:end-1)); % duration of each time segment [seconds]
if strcmp(options.method,'Planning')%assumes forecast is perfect
    date = round(864000*(date+options.Horizon/24))/864000;%%count forward by length of the horizon, rounded to nearest second
    actual_data = get_data(test_data,forecast.Timestamp,[],[]);
elseif strcmp(options.method,'Dispatch') || strcmp(options.method,'Control')   
    date = round(864000*(date+options.Resolution/24))/864000;%% count forward 1 step, rounded to nearest second
    actual_data = get_data(test_data,date(1),[],[]);
end
if strcmp(options.method,'Planning')%assumes forecast is perfect
    
    for i = 1:1:n_b
        [Tzone,Twall] = building_simulate(building(i),forecast.Weather.DrybulbC,...
            forecast.Weather.RHum,dt,forecast.Building.InternalGains(:,i),...
            forecast.Building.ExternalGains(:,i),solution.Buildings.Cooling(:,i),...
            solution.Buildings.Heating(:,i),forecast.Building.AirFlow(:,i),...
            forecast.Building.Damper(:,i),building(i).Tzone,building(i).Twall);
        building(i).Tzone = Tzone(end);
        building(i).Twall = Twall(end);
        building(i).Timestamp = date(1);
        design.Buildings(si+1:si+n_s,i) = Tzone(2:end);
    end
    for i = 1:1:n_fl
        fluid_loop(i).fluid_temperature = fluid_loop_step(gen,subnet.CoolingWater.Equipment{i},fluid_loop(i),fluid_loop(i).fluid_temperature,solution.Dispatch(2:end,:),dt);
    end
    for i = 1:1:n_g
        gen(i).CurrentState(1) = solution.Dispatch(end,i);
        gen(i).Status = solution.Dispatch(end,i)>0;
        if strcmp(gen(i).Type,{'Hydro Storage'})
            n = gen(i).QPform.Hydro.subnetNode;%dam #
            gen(i).CurrentState(2) = solution.hydroSOC(end,n);
        end
    end
    
    design.Timestamp(si+1:si+n_s) = actual_data.Timestamp;
    design.GeneratorState(si+1:si+n_s,:) = solution.Dispatch(2:end,:);
    design.LineFlows(si+1:si+n_s,:) = solution.LineFlows;
    
    %%copy select fields of actual data into saved design results
    F = fieldnames(actual_data);
    F = F(~strcmp('Timestamp',F));
    for j = 1:1:length(F)
        if isstruct(actual_data.(F{j}))
            S = fieldnames(actual_data.(F{j}));
            for i = 1:1:length(S)
                if ~isempty(actual_data.(F{j}).(S{i})) && isnumeric(actual_data.(F{j}).(S{i}))
                    design.(F{j}).(S{i})(si+1:si+n_s,:) = actual_data.(F{j}).(S{i});
                end
            end
        elseif ~isempty(actual_data.(F{j})) && isnumeric(actual_data.(F{j}))
            design.(F{j})(si+1:si+n_s,:) = actual_data.(F{j});
        end
    end
    
    if isfield(solution,'LBRelax')
        design.LBRelax(si+1:si+n_s) = solution.LBRelax;
    end
    if isfield(solution,'hydroSOC') && ~isempty(solution.hydroSOC)
        design.Hydro.SOC(si+1:si+n_s,:) = solution.hydroSOC;
        design.Hydro.spillFlow(si+1:si+n_s,:) = solution.hydroSpillFlow;
        for n = 1:1:length(subnet.Hydro.nodes)
            design.Hydro.OutFlow(si+1:si+n_s,n) = solution.LineFlows(:,subnet.Hydro.lineNumber(n));
        end
    end
    si = si + n_s; %if 24 hours, take 24 steps
elseif strcmp(options.method,'Dispatch') || strcmp(options.method,'Control')   
    si = si+1;
    
    for i = 1:1:n_g
        gen(i).CurrentState(1) = solution.Dispatch(2,i);
        gen(i).Status = solution.Dispatch(2,i)>0;
        if strcmp(gen(i).Type,'Hydro Storage')
            n = gen(i).QPform.Hydro.subnetNode;%dam #
            gen(i).CurrentState(2) = solution.hydroSOC(1,n);
        end
    end
    
    dispatch.Timestamp(si) = date(1);
    dispatch.GeneratorState(si,:) = solution.Dispatch(2,:);
    dispatch.LineFlows(si,:) = solution.LineFlows(2,:);
    predicted.Timestamp(:,si-1) = forecast.Timestamp;
    predicted.GenDisp(:,:,si-1) = solution.Dispatch(2:end,:);
    predicted.LineFlows(:,:,si-1) = solution.LineFlows;
    predicted.Buildings(:,:,si-1) = solution.Buildings.Temperature;
    predicted.fluid_loop(:,:,si-1) = solution.fluid_loop;
    
    %%copy select fields of actual data and forecasted data into saved dispatch and predicted
    F = fieldnames(actual_data);
    F = F(~strcmp('Timestamp',F));
    F = F(~strcmp('HistProf',F));
    for j = 1:1:length(F)
        if isstruct(actual_data.(F{j}))
            S = fieldnames(actual_data.(F{j}));
            for i = 1:1:length(S)
                dispatch.(F{j}).(S{i})(si,:) = actual_data.(F{j}).(S{i})(1,:);
                for k = 1:1:length(forecast.(F{j}).(S{i})(1,:))
                    predicted.(F{j}).(S{i})(si-1,:,k) = forecast.(F{j}).(S{i})(:,k)';
                end
            end
        else
            dispatch.(F{j})(si,:) = actual_data.(F{j})(1,:);
            for k = 1:1:length(forecast.(F{j})(1,:))
                predicted.(F{j})(si-1,:,k) = forecast.(F{j})(:,k)';
            end
        end
    end
    if isfield(solution,'LBRelax')
        predicted.LBRelax(si) = solution.LBRelax;
    end
    for i = 1:1:n_b
        [Tzone,Twall] = building_simulate(building(i),actual_data.Weather.DrybulbC(1),...
            actual_data.Weather.RHum(1),dt(1),forecast.Building.InternalGains(1,i),...
            forecast.Building.ExternalGains(1,i),solution.Buildings.Cooling(1,i),...
            solution.Buildings.Heating(1,i),forecast.Building.AirFlow(1,i),...
            forecast.Building.Damper(1,i),building(i).Tzone,building(i).Twall);
        building(i).Tzone = Tzone(2);
        building(i).Twall = Twall(2);
        building(i).Timestamp = date(1);
        dispatch.Buildings(si,i) = Tzone(2);
    end
    for i = 1:1:n_fl
        fluid_loop(i).fluid_temperature = fluid_loop_step(gen,subnet.CoolingWater.Equipment{i},fluid_loop(i),fluid_loop(i).fluid_temperature,solution.Dispatch(2,:),dt(1));
        dispatch.fluid_loop(si,i) = fluid_loop(i).fluid_temperature;
    end
    if isfield(solution,'hydroSOC') && ~isempty(solution.hydroSOC)
        dispatch.Hydro.SOC(si,:) = solution.hydroSOC(1,:);
        dispatch.Hydro.spillFlow(si,:) = solution.hydroSpillFlow(1,:);
        predicted.hydroSOC(:,:,si-1) = solution.hydroSOC;
        for n = 1:1:length(subnet.Hydro.nodes)
            dispatch.Hydro.OutFlow(si,n) = solution.LineFlows(1,subnet.Hydro.lineNumber(n));
        end
    end
end
if strcmp(options.method,'Control')
    %count forward in time in the control loop
        %% Real-time control 
    if any(strcmp(options.mode,{'observer';'controller';}))
        if isempty(timerfindall)
            Timers(options)
        else %do nothing;
        end
    else
        %% Virtual Plant
        D = date(1);
        while date(1)<(D+options.Resolution/24)
            OnlineLoop
        end
    end
end
end%Ends function dispatch_record
