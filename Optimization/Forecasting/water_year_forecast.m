function [wy_forecast, gen] = water_year_forecast(gen,buildings,cool_tower,subnet,options,date,wy_forecast,test_data)
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
    D = datevec(date_now);
    if D(2)<10
        year = D(1)-1;
    else
        year = D(1);
    end
    date = datenum([year 10 1 1 0 0])+[0;build_time_vector(options)/24];
    dt = (date(2:end) - date(1:end-1))*24;
    op_mat_a = load_matrices(gen,buildings,cool_tower,subnet,options,'A',dt); %build quadratic programming matrices for FitA
    op_mat_b = load_matrices(gen,buildings,cool_tower,subnet,options,'B',dt);%build quadratic programming matrices for FitB
    one_step = load_matrices(gen,buildings,cool_tower,subnet,options,'B',[]);%build quadratic programming matrices for single time step
    
    freq = 364; %period of repetition (364 = 1 year)
    res = 7;
    n_o = round(freq/res)+1;
    prev_data = get_data(test_data,linspace((date(1) - res - freq),date(1)-res,n_o)',[],[]);
    now_data = get_data(test_data,date(1),[],[]);
    future_data = get_data(test_data,date(2:end),[],[]);
    if isempty(wy_forecast)
        hydro_soc_init = zeros(1,length(subnet.Hydro.nodes));
        for n=1:1:length(subnet.Hydro.nodes)
            i = subnet.Hydro.Equipment{n};
            if strcmp(gen(i).Type,'Hydro Storage')
                gen(i).CurrentState(2) = gen(i).VariableStruct.StartWYstate*gen(i).QPform.Stor.UsableSize;
                hydro_soc_init(n) = gen(i).CurrentState(2);
            end 
        end 
        [data_t0,gen,~] = update_forecast(gen,buildings,cool_tower,subnet,options,date(1),test_data.HistProf,prev_data,now_data,future_data);
        [gen,cool_tower] = automatic_ic(gen,buildings,cool_tower,subnet,date(1),one_step,options,data_t0);% set the initial conditions
    end
    [forecast,gen,buildings] = update_forecast(gen,buildings,cool_tower,subnet,options,date(2:end),test_data.HistProf,prev_data,now_data,future_data);
    [solution,~] = dispatch_loop(gen,buildings,cool_tower,subnet,op_mat_a,op_mat_b,one_step,options,date,forecast,wy_forecast);
    solution.Timestamp = date;
    solution.hydroSOC = [hydro_soc_init;solution.hydroSOC];
    disp(strcat('Water Year Forecast Completed for ',num2str(year),':',num2str(year+1)))
%     %% Sanity Check
%     n_d = length(test_data.Hydro.OutFlow(1,:));
%     for i = 1:1:n_d
%         figure(i)
%         plot(solution.Timestamp(2:end),solution.LineFlows(:,subnet.Hydro.lineNumber(i)),'b')
%         hold on
%         plot(test_data.Timestamp(1:8737),test_data.Hydro.OutFlow(1:8737,i),'r')
%     end
%     %%---%%
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
end%Ends function water_year_forecast
