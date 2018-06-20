function plot_project(gen,buildings,fluid_loop,n_plot,dispatch,forecast,test_data,handles,mode)
if ~isempty(buildings)
    n_plot = n_plot + 1;
end
  
if strcmp(mode,'Result') %plot dispatch results in EDC
    plot_gen(gen,buildings,fluid_loop,dispatch,forecast,handles,n_plot,mode)
elseif strcmp(mode,'Dispatch') %plot dispatch results in EPC
    plot_disp(gen,buildings,fluid_loop,dispatch,forecast,handles,mode)
elseif strcmp(mode,'Building') %plot dispatch results for building temperature profile in EPC
    
    
elseif strcmp(mode,'Costs') %plot cost results in EPTother 
    plot_costs(test_data,handles)
elseif strcmp(mode,'Forecast') %plot demands on the forecasting tab
    plot_forecasts(forecast,test_data,buildings,n_plot,handles);
elseif strcmp(mode,'Demand')
    [m,n] = size(test_data.data);
    if n>length(test_data.color)
        color_vec = get(handles(3),'UserData');
        colors_plot = interp1(linspace(0,1,length(color_vec)),color_vec,linspace(0,1,n));
    else
        colors_plot = test_data.color;
    end
    set(handles(1),'ColorOrder',colors_plot)%'defaultAxesColorOrder'
    set(handles(1),'YLimMode','auto','YTickLabelMode','auto','YTickMode','auto');
    plot_demands(handles(1),test_data.Timestamp,test_data.data,[test_data.y_lab, test_data.units])
    plot_histogram(handles(2),test_data.data,[test_data.y_lab, test_data.units])   
end

end%ends function plot_project

function plot_costs(data,handles)
cla(handles(1));
cla(handles(2));
n_sys = length(data.names); %number of systems to plot
labelMonths = {};
m = {'Jan';'Feb';'Mar';'Apr';'May';'Jun';'Jul';'Aug';'Sep';'Oct';'Nov';'Dec';};
D = datevec(data.timestamp(1));
n = 0;
while datenum([D(1) D(2)+n 1])<data.timestamp(end)
    labelMonths(end+1) = m(rem(D(2)+n-1,12)+1);
    n = n+1;
end
Xpos = zeros(n,n_sys);
if n_sys == 1
    Xpos(:,1) = (.5:1:n)';
else
    for i = 1:1:n_sys
        Xpos(:,i) = ((1/n_sys^2 + (i-1)/n_sys):1:n)';
    end
end
colormap(handles(1),'summer')
for i = 1:1:n_sys
    %savebuilding color and figure out spacing/width
    bar(handles(1),Xpos(:,i),data.costs(:,:,i),'stacked','BarWidth',0.8/n_sys)
    params = set(handles(1));
    if isfield(params,'ColorOrderIndex')
        set(handles(1),'ColorOrderIndex',1);
    end
end
legend(handles(1),{'Financing Charges';'O & M Charges';'Re-start Charges';'Demand Charges';'Electric Use Charges';'Fuel Charges';})
xlim(handles(1),[0,n])
set(handles(1),'YLimMode','auto','YTickLabelMode','auto','YTickMode','auto');
ylabel(handles(1),'Cost ($)','Color','k','FontSize',14)
set(handles(1),'XTick',mean(Xpos,2),'XTickLabel', labelMonths)
colormap(handles(2),'summer')
bar(handles(2),Xpos(1,:),data.npc)
xlim(handles(2),[0,1])
ylabel(handles(2),'20 Year Net Present Cost ($)','Color','k','FontSize',14)
set(handles(2),'XTick',Xpos(1,:),'XTickLabel',data.names) 
if length(data.names)>1
    x_lab = strcat('Order is : ',data.names{1});
    for i = 2:1:n_sys
        x_lab = strcat(x_lab,' , ',data.names{i});
    end
    xlabel(handles(1),x_lab)
end
end%Ends function plot_costs

function plot_forecasts(forecast,history,buildings,n_plot,handles)
set(handles.sliderZoom,'Visible','on')
set(handles.textHour,'Visible','on'); set(handles.textWeek,'Visible','on'); set(handles.textHorizon,'Visible','on');
for i = 1:1:n_plot
    set(handles.(strcat('ForecastPlot',num2str(i))),'Visible','on');%plotting axes (y-axis on left)
    set(handles.(strcat('ForecastPlot',num2str(i),'b')),'Visible','on');%y-axis on the right
    set(handles.(strcat('ForecastName',num2str(i))),'Visible','on');%button with name
end
%% need to add heating and cooling
n_b = length(buildings);
if n_b>0 
    forecast.Demand.T = zeros(length(forecast.Building(1).E0),n_b);
    if ~isfield(forecast,'Demand') || ~isfield(forecast.Demand,'E')
        forecast.Demand.E = 0;
    end
    if ~isfield(forecast.Demand,'H')
        forecast.Demand.H = 0;
    end
    if ~isfield(forecast.Demand,'C')
        forecast.Demand.C = 0;
    end
    for i = 1:1:n_b
        forecast.Demand.E = forecast.Demand.E + forecast.Building(i).E0;
        forecast.Demand.H = forecast.Demand.H + forecast.Building(i).H0;
        forecast.Demand.C = forecast.Demand.C + forecast.Building(i).C0;
        forecast.Demand.T(:,i) = forecast.Building(i).Tzone;
    end
end

if ~isempty(history)
    history_time = (history.Timestamp - floor(history.Timestamp(1)))*24;
    forecast_time = (forecast.Timestamp - floor(history.Timestamp(1)))*24;
    timestamps = [history.Timestamp(1:nnz(history.Timestamp<forecast.Timestamp(1)));forecast.Timestamp];
else
    forecast_time = (forecast.Timestamp - floor(forecast.Timestamp(1)))*24;
    timestamps = forecast.Timestamp;
end
date_text = fixed_date_text(timestamps);
timestamps = (timestamps - floor(timestamps(1)))*24;
%% Do actual Plotting
for q = 1:1:n_plot
    h = handles.(strcat('ForecastPlot',num2str(q)));
    cla(h);
    if q==1
        text_size = 12;
    else
        text_size = 9;
    end
    y_min = 0;
    s = get(handles.(strcat('ForecastName',num2str(q))),'String');
    if strcmp(s,'Electrical')
        s = 'E';
        y_lab = 'Electrical Demand (kW)';
    elseif strcmp(s,'DistrictCurrent')
        s = 'DC';
        y_lab = 'Electrical Demand (kW)';
    elseif strcmp(s,'DistrictHeat')
        s = 'H';
        y_lab = 'Heating Demand (kW)';
    elseif strcmp(s,'DistrictCool')
        s = 'C';
        y_lab = 'Cooling Deamnd (kW)';
    elseif strcmp(s,'Hydro')
        s = 'W';
        y_lab = 'Withdrawls (1000cfs)';
    elseif strcmp(s,'Steam')
        s = 'S';
        y_lab = 'Steam Demand (1000 lbs/hr)';
    elseif strcmp(s,'BuildingTemp')
        s = 'T';
        y_lab = 'Building Temperature Setpoint (C)';
    end
    if isfield(forecast.Demand,s)
        x_tick_points = (ceil(timestamps(1)):round((timestamps(end)-timestamps(1))/12):timestamps(end));
        x_tick_label = mod(x_tick_points,24);
        x_tick_label([false,x_tick_label(2:end)==0]) = 24;

        oom = max(1,log10(max(sum(forecast.Demand.(s),2))));
        if (oom-floor(oom))==0 %count in increments of 1, 10, 100 or 1000 etc
            y_space = 10^(oom-1);
            y_max = 10^oom;
        elseif (oom-floor(oom))> 0.6990 %count in increments of 1, 10, 100 or 1000 etc
            y_space = 10^floor(oom);
            y_max = 10^ceil(oom);
        elseif (oom-floor(oom))> 0.30103 %count in increments of 5, 50, 500 or 5000 etc
            y_space = .5*10^floor(oom);
            y_max = .5*10^ceil(oom);
        else  %count in increments of 2, 20, 200 or 2000 etc
            y_space = .2*10^floor(oom);
            y_max = .2*10^ceil(oom);
        end
        if strcmp(s,'T')
            y_min = 55;
            y_max = 85;
            y_space = 5;
            forecast.Demand.T = forecast.Demand.T*9/5+32;
        end
        plot(h,forecast_time,forecast.Demand.(s),'k')
        if ~isempty(history) && isfield(history,'Demand') && isfield(history.Demand,s)
            plot(h,history_time,history.Demand.(s),'b');
        end
        xlabel(h,date_text,'Color','k','FontSize',text_size)
        ylabel(h,y_lab,'Color','k','FontSize',text_size)
        set(h,'XTick',x_tick_points,'XTickLabel', {x_tick_label})
        xlim(h,[timestamps(1), timestamps(end)])
        ylim(h,[y_min,y_max])
        set(h,'YTick',y_min:y_space:y_max,'FontSize',text_size-2)
        h2 = handles.(strcat('ForecastPlot',num2str(q),'b'));
        set(h2,'xtick',[],'xticklabel',[],'YTick',[],'YTickLabel', [])
        ylabel(h2,[]);
    end
end
end%ends function plot_forecast

function plot_gen(gen,buildings,fluid_loop,dispatch,predicted,handles,n_plot,mode)
n_g = length(gen);
n_b = length(buildings);
n_fl = length(fluid_loop);
names = cell(n_g+n_b,1);
for i = 1:1:n_g
    names(i) = {gen(i).Name};
end
for i = 1:1:n_b
    names(n_g+i) = {buildings(i).Name};
end
for i = 1:1:n_fl
    names(n_g+n_b+i) = {fluid_loop(i).location};
end
nnList = {'Electrical';'DistrictHeat';'DistrictCool';'Hydro';'DirectCurrent';'CoolingWater';'Transmission1';'Hydrogen';'LiqHydrogen';'Heating2';'BuildingTemp'};
nnAbrev = {'E';'H';'C';'W';'DC';'CW';'E1';'Hy';'LH2';'H2';'BT'};

LabelYaxis = {'Electrical Generation ';'Heating Generation ';'Cooling Generation ';'Withdrawls ';'DC Electrical Generation ';'CoolingWater Temperature ';'230kV Transmission';'Hydrogen ';'Liquid H2 ';'Steam ';'Temperature'};

for i = 1:1:n_plot
    if strcmp(get(handles.(strcat(mode,'Plot',num2str(i))),'Visible'),'on')
        han = [handles.figure1;handles.(strcat(mode,'Plot',num2str(i)));handles.(strcat(mode,'Plot',num2str(i),'b'));handles.(strcat('ResultStorageSelect',num2str(1)));];
        cla(han(2));
        cla(han(3));
        net_name = get(handles.(strcat(mode,'Name',num2str(i))),'String');
        index = nonzeros((1:11)'.*strcmp(net_name,nnList));
        ylab = LabelYaxis{index};
        net_abrv = nnAbrev{index};
        [data,data_2,timestamps,start_index] = collect_data(gen,buildings,fluid_loop,dispatch,predicted,net_abrv);
        date_text = scroll_date_text(i,timestamps,start_index);
        line = false;
        if strcmp(net_abrv,'BT') || (any(any(data_2)) && (get(handles.StackedGraph,'Value')==0 || strcmp(get(handles.uipanelMain3,'visible'),'on')))
            line = true;
        end
        if i == 1
            plot_to_gui(data,data_2,timestamps,names,date_text,net_abrv,han,line,12,1,ylab,mode)
            if start_index>1
                t_now = (timestamps(start_index) - floor(timestamps(1)))*24;
                plot(han(2),[t_now,t_now],[-1e8,1e8],'c--')  
            end
        else
            plot_to_gui(data,data_2,timestamps,names,date_text,net_abrv,han,line,9,start_index,ylab,mode)
        end
    end
end

end %ends plot_gen

function plot_disp(gen,buildings,fluid_loop,dispatch,predicted,handles,mode)
n_g = length(gen);
n_b = length(buildings);
names = cell(n_g+n_b,1);
for i = 1:1:n_g
    names(i) = {gen(i).Name};
end
for i = 1:1:n_b
    names(n_g+i) = {buildings(i).Name};
end
nnList = {'Electrical';'DistrictHeat';'DistrictCool';'Hydro';'DirectCurrent';'CoolingWater';'Transmission1';'Hydrogen';'LiqHydrogen';'Heating2';'BuildingTemp'};
nnAbrev = {'E';'H';'C';'W';'DC';'CW';'E1';'Hy';'LH2';'H2';'BT'};

LabelYaxis = {'Electrical Generation ';'Heating Generation ';'Cooling Generation ';'Withdrawls ';'DC Electrical Generation ';'CoolingWater Temperature ';'230kV Transmission';'Hydrogen ';'Liquid H2 ';'Steam ';'Temperature'};
han = [handles.figure1;handles.axesMain;handles.axesMainR;];
cla(han(2));
cla(han(3));
net_names = get(handles.popupmenuAxes,'String');
net = net_names{get(handles.popupmenuAxes,'Value')};
index = nonzeros((1:11)'.*strcmp(net,nnList));
if isempty(index)%building
   index = 11; 
end
ylab = LabelYaxis{index};
net_abrv = nnAbrev{index};
[data,data_2,timestamps,start_index] = collect_data(gen,buildings,fluid_loop,dispatch,predicted,net_abrv);
date_text = scroll_date_text(2,timestamps,start_index);
line = false;
if strcmp(net_abrv,'BT')
    line = true;
end
plot_to_gui(data,data_2,timestamps,names,date_text,net_abrv,han,line,12,start_index,ylab,mode)

end %ends plot_disp

function [data,data_2,timestamps,s0] = collect_data(gen,buildings,fluid_loop,dispatch,predicted,S)
n_g = length(gen);
n_b = length(buildings);
if ~isempty(dispatch)
    timestamps = [dispatch.Timestamp;predicted.Timestamp(2:end)];
    s0 = length(dispatch.Timestamp);
else
    timestamps = predicted.Timestamp;
    s0 = 1;
end
n_s = length(timestamps);
dt = (timestamps(2:end) - timestamps(1:end-1))*24;
data = zeros(n_s,n_g+n_b);
data_2 = zeros(n_s,n_g+n_b);
for i = 1:1:n_g
    if isfield(gen(i).QPform.output,S) 
        if ~isempty(dispatch)
            data(:,i) = [dispatch.Dispatch(:,i);predicted.Dispatch(2:end,i)];
        else
            data(:,i) = predicted.Dispatch(:,i);
        end
        switch gen(i).Type
            case 'Utility'
                %do nothing sellback acounted for in sort_solution
            case 'AC_DC'
                data(:,i) = data(:,i)*gen(i).QPform.output.(S)(1); 
            case {'Electric Generator';'Heater';'Solar';'Wind';}
                data(:,i) = data(:,i)*gen(i).QPform.output.(S);                
            case 'CHP Generator'
                %overwrite electric power with heat recovery
                if strcmp(S,'H')
                    data(:,i) = chp_heat(gen(i).QPform,data(:,i));
                end
            case 'Chiller'
                %overwrite cooling power with input
                if strcmp(S,'E') || strcmp(S,'H')
                    data(:,i) = -chill_input(gen(i),data(:,i));
                elseif strcmp(S,'CW') 
                    data(:,i) = data(:,i) + chill_input(gen(i),data(:,i));
                end
            case 'Hydrogen Generator'
                
            case 'Electrolyzer'
                
            case 'Cooling Tower'
                if strcmp(S,'E')
                    data(:,i) = -cool_tower_input(gen(i),data(:,i));
                elseif strcmp(S,'CW') 
                    data(:,i) = data(:,i)*gen(i).QPform.output.(S);  
                end
            case {'Thermal Storage';'Hydrogen Storage';'Electric Storage';}
                %compute state of charge and overwrite with power
                soc = data(:,i);
                if isfield(gen(i).VariableStruct,'MaxDOD')
                    soc = soc+ gen(i).Size*(1-gen(i).VariableStruct.MaxDOD/100); %add the unusable charge
                end
                data(2:end,i) = (soc(1:end-1) - soc(2:end))./dt; %%need to factor in efficiency
                data_2(:,i) = soc;
            case {'Hydro Storage';}
                if strcmp(S,'E')
                    data(:,i) = data(:,i)*gen(i).QPform.output.(S)(1);
                else
                    n = gen(i).QPform.Hydro.down_river;
                    flow_now = gen(i).CurrentState(1)*gen(i).QPform.Stor.Power2Flow;
                    if isempty(dispatch)
                        data(:,i) = [flow_now;predicted.LineFlows(:,n)];%Outflow is after Power Produced
                        data_2(:,i) = predicted.Dispatch(:,i); %state-of-charge
                    else
                        data(:,i) = [dispatch.LineFlows(:,n);predicted.LineFlows(:,n)];%Outflow is after Power Produced
                        data_2(:,i) = [dispatch.Dispatch(:,i);predicted.Dispatch(2:end,i)]; %state-of-charge
                    end
                end
        end
    end
end
n_b = length(buildings);
if strcmp(S,'BT')
    for i = 1:1:n_b
        if ~isempty(dispatch)
            data(:,n_g+i) = [dispatch.Buildings;predicted.Buildings.Temperature]*9/5+32;
        else
            data(:,n_g+i) = [buildings(i).Tzone;predicted.Buildings.Temperature]*9/5+32;
        end
    end
end
n_fl = length(fluid_loop);
if strcmp(S,'CW')
    for i = 1:1:n_fl
        if ~isempty(dispatch)
            data(:,n_g+i) = [dispatch.fluid_loop;predicted.fluid_loop]*9/5+32;
        else
            data(:,n_g+i) = [fluid_loop(i).fluid_temperature;predicted.fluid_loop]*9/5+32;
        end
    end
end
end%ends function collect_data

function plot_to_gui(data,data_2,timestamps,names,date_text,S,handles,line,text_size,index_start,ylab,mode)
color_vec = get(handles(1),'UserData');
colors_plot = interp1(linspace(0,1,length(color_vec)),color_vec,linspace(0,1,length(names)));
cla(handles(2));
[data,units,Ymax,Ymin,Yspace] = y_axis_scaling(data,S,index_start,line);
data_plot(handles(2),line,data,timestamps,index_start,names,colors_plot,date_text,ylab,units,text_size,Ymax,Ymin,Yspace);
plot_storage(data_2,handles,timestamps,index_start,names,ylab,colors_plot,text_size,mode)
end%Ends function plot_to_gui

function plot_storage(data_2,handles,timestamps,index_start,names,y_lab,colors_plot,text_size,mode)
%% Plot storage state-of-charge on same plot
n = length(data_2(1,:));
timestamps = (timestamps - floor(timestamps(1)))*24;
if strcmp(mode,'Result') && nnz(sum(data_2(index_start:end,:),1))>2
    if strcmp(get(handles(4),'Visible'),'off')
        stor_names = {};
        for i = 1:1:n
            if any(data_2(index_start:end,i)>0)
                stor_names(end+1) = names(i);
            end
        end
        set(handles(4),'Visible','on','value',1,'string',stor_names)
    end
    stor_num = get(handles(4),'value');%select only 1 storage at a time
    j = 0;
    for i = 1:1:n
        if any(data_2(index_start:end,i)>0)
            j = j+1;
            if j ~=stor_num
                data_2(:,i) = 0;
            end
        end
    end
elseif strcmp(mode,'Result')
    if strcmp(get(handles(4),'Visible'),'on')
        set(handles(4),'Visible','off')
    end
end
stor_names = {};
for i = 1:1:n
    if any(data_2(index_start:end,i)>0)
        stor_names(end+1) = names(i);
        L = plot(handles(3),timestamps(index_start:end),data_2(index_start:end,i));
        set(L,'Color',colors_plot(i,:),'LineStyle','-','LineWidth',2,'Marker','x','MarkerEdgeColor','k','MarkerSize',5)
    end
end

if ~isempty(stor_names)
    xlim(handles(3),get(handles(2),'Xlim'));
    set(handles(3),'xtick',[],'xticklabel',[])
    ticks = get(handles(2),'YTick');
    oom_stor = log10(max(max(data_2(index_start:end,:))));
    if (oom_stor-floor(oom_stor))==0 %count in increments of 1, 10, 100 or 1000 etc
        Ymax = 10^oom_stor;
    elseif (oom_stor-floor(oom_stor))> 0.6990 %count in increments of 1, 10, 100 or 1000 etc
        Ymax = 10^ceil(oom_stor);
    elseif (oom_stor-floor(oom_stor))> 0.30103 %count in increments of 5, 50, 500 or 5000 etc
        Ymax = .5*10^ceil(oom_stor);
    else  %count in increments of 2, 20, 200 or 2000 etc
        Ymax = .2*10^ceil(oom_stor);
    end
    pTicks = nnz(get(handles(2),'YTick')>0); % # of positive tick marks
    Yspace = Ymax/pTicks;
    negTicks = min(ticks)/(ticks(end)-ticks(end-1));
    Ymin = Yspace*negTicks;  
    ylim(handles(3),[Ymin,Ymax])
    axTickY = Ymin:Yspace:Ymax;
    set(handles(3),'YTick',axTickY,'YTickLabel', {axTickY},'FontSize',text_size-2)
    if strcmp(y_lab,'Withdrawls ')
        ylabel(handles(3),'State of Charge (1000 acre-ft)','Color','k','FontSize',text_size)
    else
        ylabel(handles(3),'State of Charge (kWh)','Color','k','FontSize',text_size)
    end
    legend(handles(3),stor_names,'Fontsize',text_size-1,'Orientation','Horizontal','Location','North','Box','off','color','none');%,'Xcolor',[1 1 1],'Ycolor',[1 1 1])
else
    set(handles(3),'xtick',[],'xticklabel',[],'YTick',[],'YTickLabel',[])
end
end%ends function plot_storage

function [data,units,Ymax,Ymin,Yspace] = y_axis_scaling(data,S,s0,line)
if ~line && any(any(data<0))
    pos = sum(max(data(s0:end-1,:),0),2);
    neg = sum(min(0,data(s0:end-1,:)),2);
    oom = max(1,.2+log10(max(pos)-min(neg)));
else
    oom = max(1,.2+log10(max(sum(data(s0:end,:),2))));
end

if oom>6.30103
    if strcmp(S,'W') == 4
        units = '1e9 cfs';
    else
        units = 'GW';
    end
    oom = oom-6;
    data = data/1e6;
elseif oom>3.30103
    if strcmp(S,'W')
        units = '1e6 cfs';
    else
        units = 'MW';
    end
    oom = oom-3;
    data = data/1000;
else
    if strcmp(S,'W')
        units = '1000 cfs';
    elseif strcmp(S,'CW') || strcmp(S,'BT')
        units = 'F';
    else
        units = 'kW';
    end
end
if strcmp(S,'BT')
    Ymax = 78;
    Ymin = 56;
    Yspace = 2;
else
    if (oom-floor(oom))==0 %count in increments of 1, 10, 100 or 1000 etc
        Yspace = 10^(oom-1);
        Ymax = 10^oom;
    elseif (oom-floor(oom))> 0.6990 
        Yspace = 10^floor(oom);
        Ymax = 10^ceil(oom);
    elseif (oom-floor(oom))> 0.30103 %count in increments of 5, 50, 500 or 5000 etc
        Yspace = .5*10^floor(oom);
        Ymax = .5*10^ceil(oom);
    else  %count in increments of 2, 20, 200 or 2000 etc
        Yspace = .2*10^floor(oom);
        Ymax = .2*10^ceil(oom);
    end

    negTicks = floor(min(min(min(0,data)))/Yspace);
    if abs(negTicks)>3
        Yspace = 2*Yspace;
        negTicks = floor(min(min(min(0,data)))/Yspace);
    end
    Ymin = Yspace*negTicks;
    if isempty(Ymin)
        Ymin = 0;
    end
end
end%Ends function Y-axis scaling

function data_plot(h,line,data,timestamps,index_start,names,colors_plot,date_text,ylab,units,text_size,Ymax,Ymin,Yspace)
n = length(data(1,:));
timestamps = (timestamps - floor(timestamps(1)))*24;
dt = timestamps(2:end) - timestamps(1:end-1);
plot_time = zeros(2*length(timestamps(index_start:end))-2,1);
plot_time(1:2:2*length(timestamps(index_start:end))-2) = [timestamps(index_start);timestamps(index_start+2:end)-.9999*dt(index_start+1:end)];
plot_time(2:2:2*length(timestamps(index_start:end))-2) = timestamps(index_start+1:end);

n_s_2 = length(data(index_start:end,1))*2-2;
step_data = zeros(n_s_2,nnz(any(data(index_start:end,:))));
if ~isempty(step_data)
    j = 0;
    name = {};
    neg_bar_names = {};
    for i = 1:1:n
        if any(data(index_start:end,i)~=0)
            j = j+1;
            step_data(1:2:end,j) = data(index_start+1:end,i);
            step_data(2:2:end,j) = data(index_start+1:end,i);
            name(end+1) = names(i);
            if any(step_data(:,j)<0)
                neg_bar_names(end+1) = names(i);
            end
        end
    end

    if line
        str2 = 'Color';
        h1 = plot(h,plot_time,step_data,'LineWidth',3);
    else
        pos_bars = max(0,step_data);
        neg_bars = min(0,step_data);
        index_keep = any(neg_bars<0);
        neg_bars = neg_bars(:,index_keep);
        str2 = 'FaceColor';
        h1 = area(h,plot_time,pos_bars,'Linestyle','none');
        if ~isempty(neg_bars)
            L = area(h,plot_time,neg_bars,'Linestyle','none');
            for c = 1:1:length(L)
                set(L(c),str2,colors_plot(strcmp(neg_bar_names(c),names),:));
            end
        end
    end
    for c = 1:1:length(h1)
        set(h1(c),str2,colors_plot(strcmp(name(c),names),:));
    end
    axTick = (ceil(timestamps(1)):round((timestamps(end)-timestamps(1))/12):timestamps(end));
    axIndex = mod(axTick,24);
    axIndex([false,axIndex(2:end)==0]) = 24;
    set(h,'XTick',axTick,'XTickLabel', {axIndex})

    xlim(h,[timestamps(index_start), timestamps(end)])
    ylim(h,[Ymin,Ymax])
    axTickY = Ymin:Yspace:Ymax;
    set(h,'YTick',axTickY,'YTickLabel', {axTickY},'FontSize',text_size-2)
    xlabel(h,date_text,'Color','k','FontSize',text_size)
    ylabel(h,strcat(ylab,{'  ('},units,')'),'Color','k','FontSize',text_size)
end
end%Ends function data_plot


function date_text = fixed_date_text(time)
months = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
day_1 = datevec(time(1));
day_2 = datevec(time(end));
n_s = length(time);
if floor(time(1)) == floor(time(round(2/3*n_s)))
    date_text = strcat(months(day_1(2)),{' '},{num2str(day_1(3))},{'  '},{num2str(day_1(1))});
elseif floor(time(1)) == floor(time(end))-1% two days
    date_text = strcat(months(day_1(2)),{' '},{num2str(day_1(3))},{'  '},{num2str(day_1(1))},{'                       '},months(day_2(2)),{' '},{num2str(day_2(3))},{'  '},{num2str(day_2(1))});
else %many days
    date_text = strcat(months(day_1(2)),{' '},{num2str(day_1(3))},{'  '},{num2str(day_1(1))},{'   through  '},months(day_2(2)),{' '},{num2str(day_2(3))},{'  '},{num2str(day_2(1))});
end
end%Ends function fixed_date_text

function date_text = scroll_date_text(i,time,start_index)
%% Make text strings to scroll across bottom axis
months = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
vector_date = datevec(time(start_index));
day_prev = datevec(time(1)-1);
day_next = datevec(time(1)+1);
dt_1 = time(2) - time(1);
if i == 1
    date_text = ' ';
    if start_index<1/dt_1
        b = ceil(30*(1/dt_1-start_index)/(1/dt_1));
        for j = 1:1:b
            date_text = strcat(date_text,{' '});
        end
    end
    if nnz(time<floor(time(start_index)))>0.3*length(time) %more than 30% of window is previous day
        date_text = strcat(date_text,months(day_prev(2)),{' '},{num2str(day_prev(3))},{'  '},{num2str(day_prev(1))},{'                                '});
    else  %append spaces
        c = ceil(30*nnz(time<floor(time(start_index)))/length(time));
        for j = 1:1:c
            date_text = strcat(date_text,{' '});
        end
    end
    date_text = strcat(date_text,months(vector_date(2)),{' '},{num2str(vector_date(3))},{'  '},{num2str(vector_date(1))});
    if nnz(time>datenum([day_next(1),day_next(2),day_next(3)]))>0.3*length(time) %more than 30% of window is next day
        date_text = strcat(date_text,{'                         '},months(day_next(2)),{' '},{num2str(day_next(3))},{'  '},{num2str(day_next(1))});
    else  %append spaces
            date_text = strcat(date_text,{'                                      '});
    end
else
    date_text = strcat('  ',months(vector_date(2)),{' '},{num2str(vector_date(3))},{'  '},{num2str(vector_date(1))});
    if nnz(time>datenum([day_next(1),day_next(2),day_next(3)]))>0.3*length(time) %more than 30% of window is next day
        date_text = strcat(date_text,{'                      '},months(day_next(2)),{' '},{num2str(day_next(3))},{'  '},{num2str(day_next(1))});
    else
        date_text = strcat(date_text,{'                                '});
    end
end

end%ends function scroll_date_text

function plot_demands(h,timestamps,data,Ylab)
cla(h);
start_date = datevec(timestamps(1));
days = max(1,round(timestamps(end)-timestamps(1)));
monthvec = [0 31 59 90 120 151 181 212 243 273 304 334 365];
leapyear = mod((start_date(1)-2004),4)==0;%if it is a leap year it will have a remainder of 0
if leapyear
    monthvec(3:end) = monthvec(3:end)+1;
end
if days>sum(monthvec) %if you are plotting multiple years make the ticks in years
    end_date = datevec(timestamps(end));
    x_tick_points = timestamps(1)+linspace(0,round(timestamps(end)-timestamps(1)),end_date(1)-start_date(1)+1);
elseif days==sum(monthvec) %if you are plotting a year, make the ticks in months
    x_tick_points = timestamps(1)+[monthvec(start_date(2):end),monthvec(12)+monthvec(1:start_date(2))]-monthvec(start_date(2))-start_date(3);%make sure to add the right number of days given the month and day of the year.
elseif days > 31
    [y,m1,~] = datevec(timestamps(1));
    [~,m2,~] = datevec(timestamps(end));
    month=0;
    for i=1:1:(m2-m1)
        d = datenum([y,m1+i,1])-datenum([y,m1+i-1,1]);%days in month
        month(end+1) = month(end) + d;
    end
    x_tick_points = timestamps(1)+month;
elseif days>1
    x_tick_points = timestamps(1)+linspace(0,days,days+1);
else
    hours = floor(24*(timestamps(end)-timestamps(1)));
    x_tick_points = timestamps(1)+linspace(0,hours/24,floor(hours/3)+1);  
end

[m,n] = size(data);
% ch = get(h,'Children');
% for i = 1:1:length(ch)
%     if iscell(ch)
%         delete(ch{i})
%     else
%         delete(ch(i))
%     end
% end
plot(h,timestamps,data);
% for i = 2:1:n
%     plot(h,timestamps,data(:,i));
% end
ylabel(h,Ylab)
xlim(h,[timestamps(1) timestamps(end)])
set(h,'xtick',x_tick_points)
month_label = ['January  '; 'February '; 'March    '; 'April    '; 'May      '; 'June     '; 'July     '; 'August   '; 'September'; 'October  ';'November ' ;'December ';];
if days>366
    datetick(h,'x','yyyy','keeplimits')
    xlabel(h,'Year') 
elseif days==365|| days==366
    datetick(h,'x','mmmdd','keeplimits')
    xlabel(h,num2str(start_date(1)))
elseif days>=28 && days<=31
    datetick(h,'x','dd','keeplimits')
    xlabel(h,strcat(month_label(start_date(2),:),'  ',num2str(start_date(1))))
elseif days==7
    datetick(h,'x','dd','keeplimits')
    xlabel(h,strcat(month_label(start_date(2),:),'  ',num2str(start_date(1))))
elseif days ==1
    datetick(h,'x','HH','keeplimits','keepticks')
    xlabel(h,strcat(['Hours of ', month_label(start_date(2),:), num2str(start_date(3))]))
end
end %Ends function plot_demands

function plot_histogram(h,data,x_lab)
cla(h);
n = length(data);
sorted_data = sort(data);
min_val = max(1,floor(0.01*n));
max_val = ceil(0.99*n);
oom = log10(sorted_data(max_val)-sorted_data(min_val));
if (oom-floor(oom))==0 %count in increments of 1, 10, 100 or 1000 etc
    Xspace = 10^(oom-1);
elseif (oom-floor(oom))> 0.6990 %count in increments of 1, 10, 100 or 1000 etc
    Xspace = 10^floor(oom);
elseif (oom-floor(oom))> 0.30103 %count in increments of 5, 50, 500 or 5000 etc
    Xspace = .5*10^floor(oom);
else  %count in increments of 2, 20, 200 or 2000 etc
    Xspace = .2*10^floor(oom);
end
hist = [];
label = {};
min_val = ceil(sorted_data(min_val)/Xspace)*Xspace;
hist(end+1) = nnz(sorted_data<=min_val);
label(end+1) = {strcat('<',num2str(min_val))};
while min_val<sorted_data(max_val)
    min_val = min_val + Xspace;
    hist(end+1) = (nnz(sorted_data<=min_val) - sum(hist(1:end)));
    label(end+1) = {strcat(num2str(min_val-Xspace),'--',num2str(min_val))};
end
hist = hist/n*100;
bar(h,hist)
ylabel(h,'Percent of Data Within Range')
% set(h,'XTickLabel', label)
set(h,'XTickLabel','')
xlim(h,[0.5,length(hist)+.5]);
pos = get(h,'position');
t = text(0,0,x_lab);
set(t,'Parent',h,'Units','characters','HorizontalAlignment','center','Position',[pos(3)/2,-6,0]);
% Place the text labels
for i = 1:length(hist)
    t = text(0,0,label{i});
    xpos = i*pos(3)/length(hist)- 0.5*pos(3)/length(hist);
    set(t,'Parent',h,'Units','characters','HorizontalAlignment','right','VerticalAlignment','top','Rotation',45,'Position',[xpos,0,0]);
end
end %Ends function plot_histogram