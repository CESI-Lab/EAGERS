function optimize_plant_size
%OPTIMIZE_PLANT_SIZE Optimize the size of all selected generators in the system
global testSystems SYSINDEX TestData
% Find the generators in the current plant that are searchable, i.e. that
% can be considered for resizing.
handles = guihandles;
use_design_day = get(handles.DesignDay,'Value') == 1;
n_years = str2double(get(handles.NPC_Years,'String'));
range = get(handles.uitableCosts, 'Data');
i_searchable = [];
lower_bound = [];
upper_bound = [];
for i = 1:1:length(range(:,1))
    if range{i,8}
        i_searchable(end+1) = i;
        lower_bound(end+1) = range{i,9};
        upper_bound(end+1) = range{i,10};
    end
end

%%add some condition hear to go into discrete sizing mode.
%% edit the variable multiple in the call of design test to change # of generators
%% design_test(testSystems(SYSINDEX), gen_size, multiple, i_searchable, TestData, use_design_day, n_years);


if length(i_searchable) == 1
    optimize_component_size(i_searchable,lower_bound,upper_bound);
elseif length(i_searchable)>1
    % Give user feedback.
    fprintf('System optimization in progress. Please wait...\n')

    % Use particle swarm to optimize the size of each generator by minimizing
    % the net present cost of the entire plant.
    %
    % objective function:       design_test
    % optimization variable:    gen_size
    % fixed:                    curr_plant, i_searchable, TestData,
    %                           use_design_day, n_years
    % Optimization parameters.
    plot_fcn_list = {@psw_plot_best_obj_vals,@psw_plot_iter_time,@psw_plot_best_combos, @psw_plot_combos};
    options = optimoptions('particleswarm', 'Display','iter','SwarmSize', 20,'OutputFcn', @psw_output,'PlotFcn', plot_fcn_list);
    
    multiple = ones(length(i_searchable),1);
    obj_func = @(gen_size)design_test(testSystems(SYSINDEX), gen_size, multiple, i_searchable, TestData, use_design_day, n_years);
    [optim_size, ~, ~, ~] = particleswarm(obj_func, length(i_searchable), lower_bound, upper_bound, options);

    % update the plant with optimal generator sizes
    [testSystems(SYSINDEX).Generator,testSystems(SYSINDEX).Costs.Equipment,testSystems(SYSINDEX).Network] = design_resize(testSystems(SYSINDEX).Generator,testSystems(SYSINDEX).Costs.Equipment,testSystems(SYSINDEX).Network,optim_size,i_searchable);
end
MainScreen1('popupmenu_Callback')
end%ends function optimize_plant_size