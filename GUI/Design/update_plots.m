function update_plots
global testSystems TestData SYSINDEX
handles = guihandles;
ax_num = 1;
list = get(handles.popupmenuAxes,'String');
item = list{get(handles.popupmenuAxes,'Value')};
h_plot = [handles.axesMain;handles.axesCumulative;handles.figure1];
cla(handles.axesMainR);
if strcmp(item,'Monthly Costs')
    set(handles.sliderZoom1,'Visible','off');set(handles.sliderDate1,'Visible','off');set(handles.textDate1,'Visible','off');
    set(handles.textDay1,'Visible','off'); set(handles.textAllData1,'Visible','off'); set(handles.textHorizon1,'Visible','off'); 
    years = str2double(get(handles.NPC_Years,'String'));%choose years in GUI
    design_day = get(handles.DesignDay,'Value');
    base = get(handles.popupmenuBaseline,'Value');%Which project is the baseline, determined  by GUI
    [testSystems,test_data.timestamp,test_data.costs,test_data.npc,monthly_costs] = run_planning(testSystems,TestData,years,design_day,base);
    n_sys = length(testSystems);
    test_data.names = cell(n_sys,1);
    for i = 1:1:n_sys%Run through list of projects
        test_data.names(i) = {testSystems(i).Name};
        testSystems(i).Costs.Design = test_data.costs(:,:,i);
        testSystems(i).Costs.NPC = test_data.npc(:,i);
        testSystems(i).Design.Timestamp = test_data.timestamp(:,i);
        testSystems(i).Costs.ProjectedMonthlyCosts = monthly_costs(:,i);
    end
    color_vec = get(handles.figure1,'UserData');
    colors_plot = interp1(linspace(0,1,length(color_vec)),color_vec,linspace(0,1,length(test_data.costs(1,:,1))));
    set(handles.axesMain,'ColorOrder',colors_plot)%'defaultAxesColorOrder'
    plot_project([],[],[],[],[],[],test_data,h_plot,'Costs')
else
    if get(handles.sliderZoom1,'Value')== get(handles.sliderZoom1,'Max')
        set(handles.sliderDate1,'Visible','off');set(handles.textDate1,'Visible','off');
    elseif strcmp(get(handles.sliderDate1,'Visible'),'off')
        set(handles.sliderDate1,'Visible','on');set(handles.textDate1,'Visible','on');
    end
    set(handles.sliderZoom1,'Visible','on');
    set(handles.textDay1,'Visible','on'); set(handles.textAllData1,'Visible','on'); set(handles.textHorizon1,'Visible','on');
    
    if strcmp(get(handles.sliderZoom1,'Visible'),'off')
        set(handles.(strcat('sliderZoom',num2str(ax_num))),'Visible','on','value',1);set(handles.(strcat('sliderDate',num2str(ax_num))),'Visible','on','value',0);
        set(handles.(strcat('textDay',num2str(ax_num))),'Visible','on'); set(handles.(strcat('textAllData',num2str(ax_num))),'Visible','on'); set(handles.(strcat('textHorizon',num2str(ax_num))),'Visible','on');
    end
    %find the current date
    date_start = TestData.Timestamp(1) + get(handles.(strcat('sliderDate',num2str(ax_num))),'Value');
    max_slider = get(handles.(strcat('sliderZoom',num2str(ax_num))),'Max');
    if get(handles.(strcat('sliderZoom',num2str(ax_num))),'Value')==max_slider
        date_end = TestData.Timestamp(end);
    elseif floor(get(handles.(strcat('sliderZoom',num2str(ax_num))),'Value'))<1
        date_end = min(TestData.Timestamp(end),date_start + 1);
    elseif floor(get(handles.(strcat('sliderZoom',num2str(ax_num))),'Value'))<2
        date_end = min(TestData.Timestamp(end),date_start + 7);
    elseif floor(get(handles.(strcat('sliderZoom',num2str(ax_num))),'Value'))<3
        date_end = min(TestData.Timestamp(end),date_start + 31);
    elseif floor(get(handles.(strcat('sliderZoom',num2str(ax_num))),'Value'))<4
        date_end = min(TestData.Timestamp(end),date_start + 365);
    end
    switch item
        case {'Electrical Demand';'Cooling Demand';'Heating Demand';'Water Demand';'NonHVACelectric';'InternalGains';}
            x_i = nnz(TestData.Timestamp<=date_start);
            x_f = nnz(TestData.Timestamp<=date_end);
            test_data.Timestamp = TestData.Timestamp(x_i:x_f);
            test_data.y_lab = item;
            test_data.units = ' (kW)';
            test_data.color = [0 0 0];
            switch item
                case 'Electrical Demand'
                    test_data.data = TestData.Demand.E(x_i:x_f,:);
                case 'Cooling Demand'
                    test_data.data = TestData.Demand.C(x_i:x_f,:);
                    test_data.color = [0 0 1];
                case 'Heating Demand'
                    test_data.data = TestData.Demand.H(x_i:x_f,:);
                    test_data.color = [1 0 0];
                case 'Water Demand'
                    test_data.data = TestData.Demand.W(x_i:x_f,:);
                    test_data.units = ' (1000 CFS)';
                    test_data.color = [0 1 1];
                case 'NonHVACelectric'
                    test_data.data = TestData.Building.(item)(x_i:x_f,:);
                case 'InternalGains'
                    test_data.data = TestData.Building.(item)(x_i:x_f,:);
            end
            plot_project([],[],[],[],[],[],test_data,h_plot,'Demand')
        otherwise %plot  equipment schedule or building temperature profile
            %run dispatch optimization for window
            buildings = [];
            if isfield(testSystems(SYSINDEX),'Building')   
                buildings = testSystems(SYSINDEX).Building;
            end
            fluid_loop = [];
            if isfield(testSystems(SYSINDEX),'fluid_loop') 
                fluid_loop = testSystems(SYSINDEX).fluid_loop;
            end
            subnet = testSystems(SYSINDEX).Network;
            dispatch = [];op_mat_a = [];op_mat_b = [];one_step = []; predicted = [];
            if isfield(testSystems(SYSINDEX),'subnet') 
                subnet = testSystems(SYSINDEX).subnet;
                dispatch = testSystems(SYSINDEX).Dispatch;
                predicted = testSystems(SYSINDEX).Predicted;
                op_mat_a = testSystems(SYSINDEX).OpMatA;
                op_mat_b = testSystems(SYSINDEX).OpMatB;
                one_step = testSystems(SYSINDEX).OneStep;
            end
            [solution,history,gen,buildings,fluid_loop] = planning_slider_callback(...
                TestData,date_start,testSystems(SYSINDEX).Generator,buildings,fluid_loop,...
                subnet,op_mat_a,op_mat_b,one_step,dispatch,predicted,testSystems(SYSINDEX).optimoptions);
            plot_project(gen,buildings,fluid_loop,[],history,solution,[],handles,'Dispatch')
    end
end