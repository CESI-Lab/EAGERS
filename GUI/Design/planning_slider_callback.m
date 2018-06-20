% --- Executes on slider movement.
function [solution,history,gen,buildings,cool_tower] = planning_slider_callback(test_data,date,gen,buildings,cool_tower,subnet,op_mat_a,op_mat_b,one_step,dispatch,predicted,options)
handles = guihandles;
if strcmp(options.solver,'NREL')
    buildings = [];
    cool_tower = [];
    buildings.Name = 'Building';
    buildings.Tzone = 22;
%     buildings.QPform.Location.Longitude = 105; buildings.QPform.Location.Latitude = 40; buildings.QPform.Location.TimeZone = -6;
    solution = SingleOptimizationNREL(test_data.Timestamp(1) + get(handles.sliderStartDate,'Value'));
    plot_project(gen,buildings,[],3,[],solution,[],handles,'Result')
else
    if isempty(dispatch) || ~isfield(dispatch,'Timestamp') || isempty(dispatch.Timestamp) || min(abs(dispatch.Timestamp-(test_data.Timestamp(1) + get(handles.sliderStartDate,'Value'))))>=options.Resolution/24
        date = date +[0;build_time_vector(options)/24];
        test_data = update_test_data(test_data,[],gen,options);
        if isempty(dispatch)
            [gen,buildings,cool_tower,subnet,op_mat_a,op_mat_b,one_step,~] = initialize_optimization(gen,buildings,cool_tower,subnet,options,test_data);
        end
        [wy_forecast,gen] = water_year_forecast(gen,buildings,cool_tower,subnet,options,date,[],test_data);%if october 1st,Run a yearly forecast for hydrology
        freq = 1; %period of repetition (1 = 1 day)
        res = options.Resolution/24;
        n_o = round(freq/res)+1;
        prev_data = get_data(test_data,linspace((date(1) - res - freq),date(1)-res,n_o)',[],[]);
        now_data = get_data(test_data,date(1),[],[]);
        future_data = get_data(test_data,date(2:end),[],[]);
        if strcmp(options.solver,'ANN')
            options.solver = 'quadprog';
            [solution,forecast,gen,buildings,cool_tower] = single_optimization(date,[],gen,buildings,cool_tower,options,dispatch,subnet,op_mat_a,op_mat_b,one_step,test_data.HistProf,wy_forecast,prev_data,now_data,future_data);
            options.solver = 'ANN';
        else
            [solution,forecast,gen,buildings,cool_tower] = single_optimization(date,[],gen,buildings,cool_tower,options,dispatch,subnet,op_mat_a,op_mat_b,one_step,test_data.HistProf,wy_forecast,prev_data,now_data,future_data);
        end
        history = [];
    else
        date = test_data.Timestamp(1) + get(handles.sliderStartDate,'Value');
        Si = max((1:1:length(dispatch.Timestamp))'.*(dispatch.Timestamp<=date & dispatch.Timestamp>0)); %index preceeding current step
        if (dispatch.Timestamp(Si+1)>0 && (dispatch.Timestamp(Si+1)-date)<(date-dispatch.Timestamp(Si)))
            Si = Si+1; %The next time step is actually closer
        elseif predicted.Timestamp(1,Si)==0
            Si = Si - 1;
        end
        date = predicted.Timestamp(:,Si);
        solution.Dispatch = predicted.GenDisp(:,:,Si);
        solution.LineFlows = predicted.LineFlows(:,:,Si);
        solution.Buildings.Temperature = predicted.Buildings(:,:,Si);
        solution.hydroSOC = predicted.hydroSOC(:,:,Si);
        backSteps = min(Si-1,options.Horizon/options.Resolution);
        history.Dispatch = dispatch.GeneratorState(Si-backSteps:Si,:);
        history.LineFlows = dispatch.LineFlows(Si-backSteps:Si,:);
        history.Buildings = dispatch.Buildings(Si-backSteps:Si,:);
        history.hydroSOC = dispatch.hydroSOC(Si-backSteps:Si,:);
        history.Timestamp = dispatch.Timestamp(Si-backSteps:Si,:);
    end    
    solution.Timestamp = date;
end
