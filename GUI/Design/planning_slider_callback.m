% --- Executes on slider movement.
function [solution,history,gen,buildings,fluid_loop] = planning_slider_callback(test_data,date,gen,buildings,fluid_loop,subnet,op_mat_a,op_mat_b,one_step,dispatch,predicted,options)
handles = guihandles;
if strcmp(options.solver,'NREL')
    solution = single_optimization_NREL(gen,buildings,options,test_data,date);
    buildings.Tzone = solution.Buildings.Temperature(1);
    plot_project(gen,buildings,[],3,[],solution,[],handles,'Result')
else
    if isempty(dispatch) || ~isfield(dispatch,'Timestamp') || isempty(dispatch.Timestamp) || min(abs(dispatch.Timestamp-(test_data.Timestamp(1) + get(handles.sliderStartDate,'Value'))))>=options.Resolution/24
        date = date +[0;build_time_vector(options)/24];
        test_data = update_test_data(test_data,[],[],options);
        if isempty(dispatch)
            [gen,buildings,fluid_loop,subnet,op_mat_a,op_mat_b,one_step,~,~,dispatch,~,~,~] = initialize_optimization(gen,buildings,fluid_loop,subnet,options,test_data,date(1));
        end
        [wy_forecast,gen] = water_year_forecast(gen,buildings,fluid_loop,subnet,options,date,[],test_data);%if october 1st,Run a yearly forecast for hydrology

        if strcmp(options.solver,'ANN')
            options.solver = 'quadprog';
            [solution,forecast,gen,buildings,fluid_loop] = single_optimization(date,[],gen,buildings,fluid_loop,options,dispatch,subnet,op_mat_a,op_mat_b,one_step,test_data,wy_forecast);
            options.solver = 'ANN';
        else
            [solution,forecast,gen,buildings,fluid_loop] = single_optimization(date,[],gen,buildings,fluid_loop,options,dispatch,subnet,op_mat_a,op_mat_b,one_step,test_data,wy_forecast);
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
