% --- Executes on slider movement.
function slider_callback(index)
global Plant TestData
handles = guihandles;
if strcmp(Plant.optimoptions.solver,'NREL')
    Plant.Building = [];
    Plant.Building.Name = 'Building';
    Plant.Building.Tzone = 22;
%     Plant.Building.QPform.Location.Longitude = 105; Plant.Building.QPform.Location.Latitude = 40; Plant.Building.QPform.Location.TimeZone = -6;
    solution = SingleOptimizationNREL(TestData.Timestamp(1) + get(handles.sliderStartDate,'Value'));
    plot_project(Plant.Generator,Plant.Building,[],3,[],solution,[],handles,'Result')
else
    if ~isfield(Plant,'Dispatch') || isempty(Plant.Dispatch) || ~isfield(Plant.Dispatch,'Timestamp') || isempty(Plant.Dispatch.Timestamp) || min(abs(Plant.Dispatch.Timestamp-(TestData.Timestamp(1) + get(handles.sliderStartDate,'Value'))))>=Plant.optimoptions.Resolution/24
        date = TestData.Timestamp(1) + get(handles.sliderStartDate,'Value')+[0;build_time_vector(Plant.optimoptions)/24];
        if ~isfield(Plant,'Building')
            Plant.Building = [];
        end
        if ~isfield(Plant,'fluid_loop') 
            Plant.fluid_loop = [];
        end
        if isfield(Plant,'Data') 
            data = Plant.Data;
        else
            data = [];
        end
        TestData = update_test_data(TestData,data,Plant.Generator,Plant.optimoptions);
        if ~isfield(Plant,'Dispatch') || ~isempty(Plant.Dispatch)
            [Plant.Generator,Plant.Building,Plant.fluid_loop,Plant.subNet,Plant.OpMatA,Plant.OpMatB,Plant.OneStep,~] = initialize_optimization(Plant.Generator,Plant.Building,Plant.fluid_loop,Plant.Network,Plant.optimoptions,TestData);
            Plant.Dispatch = [];
        end
        [wy_forecast,Plant.Generator] = water_year_forecast(Plant.Generator,Plant.Building,Plant.fluid_loop,Plant.subNet,Plant.optimoptions,date,[],TestData);%if october 1st,Run a yearly forecast for hydrology
        dispatch = Plant.Dispatch;
        freq = 1; %period of repetition (1 = 1 day)
        res = Plant.optimoptions.Resolution/24;
        n_o = round(freq/res)+1;
        prev_data = get_data(TestData,linspace((date(1) - res - freq),date(1)-res,n_o)',[],[]);
        now_data = get_data(TestData,date(1),[],[]);
        future_data = get_data(TestData,date(2:end),[],[]);
        [prev_data,now_data] = water_dispatch_history(prev_data,now_data,dispatch,Plant.subNet);
        if strcmp(Plant.optimoptions.solver,'ANN')
            Plant.optimoptions.solver = 'quadprog';
            [solution,forecast,Plant.Generator,Plant.Building,Plant.fluid_loop] = single_optimization(date,[],Plant.Generator,Plant.Building,Plant.fluid_loop,Plant.optimoptions,dispatch,Plant.subNet,Plant.OpMatA,Plant.OpMatB,Plant.OneStep,TestData.HistProf,wy_forecast,prev_data,now_data,future_data);
            Plant.optimoptions.solver = 'ANN';
        else
            [solution,forecast,Plant.Generator,Plant.Building,Plant.fluid_loop] = single_optimization(date,[],Plant.Generator,Plant.Building,Plant.fluid_loop,Plant.optimoptions,dispatch,Plant.subNet,Plant.OpMatA,Plant.OpMatB,Plant.OneStep,TestData.HistProf,wy_forecast,prev_data,now_data,future_data);
        end
        history = [];
    else
        d = TestData.Timestamp(1) + get(handles.sliderStartDate,'Value');
        s_i = max((1:1:length(Plant.Dispatch.Timestamp))'.*(Plant.Dispatch.Timestamp<=d & Plant.Dispatch.Timestamp>0)); %index preceeding current step
        if (Plant.Dispatch.Timestamp(s_i+1)>0 && (Plant.Dispatch.Timestamp(s_i+1)-d)<(d-Plant.Dispatch.Timestamp(s_i)))
            s_i = s_i+1; %The next time step is actually closer
        elseif Plant.Predicted.Timestamp(1,s_i)==0
            s_i = s_i - 1;
        end
        date = [d;Plant.Predicted.Timestamp(:,s_i)];
        solution.Dispatch = [zeros(1,length(Plant.Predicted.GenDisp(1,:,s_i)));Plant.Predicted.GenDisp(:,:,s_i);];
        solution.LineFlows = Plant.Predicted.LineFlows(:,:,s_i);
        solution.Buildings.Temperature = Plant.Predicted.Buildings(:,:,s_i);
        solution.fluid_loop = Plant.Predicted.fluid_loop(:,:,s_i);
        solution.hydroSOC = Plant.Predicted.hydroSOC(:,:,s_i);
        backSteps = min(s_i-1,Plant.optimoptions.Horizon/Plant.optimoptions.Resolution);
        history.Dispatch = Plant.Dispatch.GeneratorState(s_i-backSteps:s_i,:);
        history.LineFlows = Plant.Dispatch.LineFlows(s_i-backSteps:s_i,:);
        history.Buildings = Plant.Dispatch.Buildings(s_i-backSteps:s_i,:);
        history.fluid_loop = Plant.Dispatch.fluid_loop(s_i-backSteps:s_i,:);
        history.hydroSOC = Plant.Dispatch.hydroSOC(s_i-backSteps:s_i,:);
        history.Timestamp = Plant.Dispatch.Timestamp(s_i-backSteps:s_i,:);
    end    
    solution.Timestamp = date;
    gen = Plant.Generator;
    mode = 'Result';
    if ~isempty(index)
        gen = gen(index);
        solution.Dispatch = solution.Dispatch(:,index);
        if ~isempty(history)
            history.Dispatch = history.Dispatch(:,index);
        end
        mode = 'Forecast';
    end
    n_plot = length(fieldnames(Plant.subNet));
    plot_project(gen,Plant.Building,Plant.fluid_loop,n_plot,history,solution,[],handles,mode)
end
