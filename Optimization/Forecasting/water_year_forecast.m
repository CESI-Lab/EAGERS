function [wy_forecast, gen] = water_year_forecast(gen,buildings,fluid_loop,subnet,options,date,wy_forecast,test_data)
hydroforecast = false;
for i = 1:1:length(gen)
    if strcmp(gen(i).Type,'Hydro Storage')
        hydroforecast = true;
    end
end
if ~isempty(wy_forecast)
    if wy_forecast.Timestamp(end)>=date(end)
        hydroforecast = false;
    end
end

if hydroforecast
    %first create the yearly dispatch data 
    % i.e. run Dispatch loop with updated information
    %these will be used for set points in the actual dispatch
    date_now = date(1);
    options.Horizon = 364*24; %Yearly Horizon
    options.Resolution = 7*24; %Week Resolution
    s_date = datevec(date(1));
    e_date = datevec(date(end));
    if s_date(2)<10 && e_date(2)>=10
        year = s_date(1)-1;
        dy1 = datenum([s_date(1),10,1,0,0,0]);
        options.Horizon = (53+ceil((date(end)-dy1)/7))*7*24;%ensure that water year forecast goes a week beyond the final date
    elseif s_date(2)<10
        year = s_date(1)-1;
    else
        year = s_date(1);
    end
    date = datenum([year 10 1 1 0 0])+[0;build_time_vector(options)/24];
    dt = (date(2:end) - date(1:end-1))*24;
    op_mat_a = load_matrices(gen,buildings,fluid_loop,subnet,options,'A',dt); %build quadratic programming matrices for FitA
    op_mat_b = load_matrices(gen,buildings,fluid_loop,subnet,options,'B',dt);%build quadratic programming matrices for FitB
    one_step = load_matrices(gen,buildings,fluid_loop,subnet,options,'B',[]);%build quadratic programming matrices for single time step
    
    if isempty(wy_forecast)
        for n=1:1:length(subnet.Hydro.nodes)
            i = subnet.Hydro.Equipment{n};
            if strcmp(gen(i).Type,'Hydro Storage')
                gen(i).CurrentState(2) = gen(i).VariableStruct.StartWYstate*gen(i).QPform.Stor.UsableSize;
            end 
        end 
        [data_t0,gen,~] = update_forecast(gen,buildings,fluid_loop,subnet,options,date(1),test_data,[]);
        [gen,fluid_loop] = automatic_ic(gen,buildings,fluid_loop,subnet,date(1),one_step,options,data_t0);% set the initial conditions       
    end
    
    if date(end)<=test_data.Timestamp(end)
        [forecast,gen,buildings] = update_forecast(gen,buildings,fluid_loop,subnet,options,date(2:end),test_data,[]);
        [solution,~] = dispatch_loop(gen,buildings,fluid_loop,subnet,op_mat_a,op_mat_b,one_step,options,date,forecast,wy_forecast);
        solution.Timestamp = date;
        hydro_soc_init = zeros(1,length(subnet.Hydro.nodes));
        for n=1:1:length(subnet.Hydro.nodes)
            i = subnet.Hydro.Equipment{n};
            if strcmp(gen(i).Type,'Hydro Storage')
                hydro_soc_init(n) = gen(i).CurrentState(2);
            end 
        end 
        solution.hydroSOC = [hydro_soc_init;solution.hydroSOC];
        disp(strcat('Water Year Forecast Completed for ',num2str(year),':',num2str(year+1)))
        if isempty(wy_forecast)
            for i = 1:1:length(gen)
                if strcmp(gen(i).Type,'Hydro Storage')
                    n = gen(i).QPform.Hydro.subnetNode;%dam #
                    gen(i).CurrentState(2) = interp1(solution.Timestamp,solution.hydroSOC(:,n),date_now);
                end
            end
        end
        wy_forecast = solution;
    end
end
end%Ends function water_year_forecast
