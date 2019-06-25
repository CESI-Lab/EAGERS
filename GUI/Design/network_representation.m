function handles = network_representation(sys,handles)
%% create diagram of equipment and network
%%collect networks, nodes and equipment
network = sys.Network;
buildings = [];
if isfield(sys,'Building') && ~isempty(sys.Building)
    buildings = sys.Building;
end
[gen, buildings] = update_qpform_all(sys.Generator,buildings,network,1);% updates the QPform field in all generators and buildings
[subnet, gen, buildings] = load_network(network,gen,buildings);

n_g = length(gen);
nodes = length(network);
gen_names = cell(n_g,1);
for i = 1:1:n_g
    gen_names(i,1) = {gen(i).Name};
end
[l_norm,node_names,h_norm] = normalized_node_location(network);
nn_list = {'Electrical';'DistrictHeat';'DistrictCool';'Hydro';'DirectCurrent';'CoolingWater';'Transmission1';'Hydrogen';'LiqHydrogen';'Heating2';};
nn_abrev = {'E';'H';'C';'W';'DC';'CW';'E1';'Hy';'LH2';'H2';};
nn_color = [.9 .9 .9; 1 .6 .6; .3 .75 .93; .85 .7 .1; .76 .87 .78; 0 1 1; 1 .6 .78; 1 .4 .4; 0 .8 1; 1 .2 .4;];
network_names = fieldnames(network);
network_names = network_names(~strcmp('name',network_names));
network_names = network_names(~strcmp('Equipment',network_names));
network_names = network_names(~strcmp('Location',network_names));
[~,reorder] = ismember(network_names,nn_list);
nn_abrev = nn_abrev(nonzeros(reorder));
nn_color = nn_color(nonzeros(reorder),:);
quote='''';
%create or update network select menu
if isfield(handles,'network_select')
    set(handles.network_select,'String',network_names,'Position', [80 46 20 2]);
else
    callback = eval(strcat('@(hObject,eventdata)MainScreen1(',quote,'network_select',quote,',hObject,eventdata,guidata(hObject))'));
    handles.network_select = uicontrol('Style', 'popupmenu', 'String',network_names, ...
        'Units','characters','FontSize', 10,'Position', [80 46 20 2],...
        'BackgroundColor','white','UserData',1,'Tag','network_select',...
        'Parent', handles.uipanelMain1,'Callback',callback,'Visible','on');%popup menu with node name and all equipment
end

%creater container for network view of each network, put a pup_up menu with
%the equipment at that node on that network. Put an axes behind it to draw
%the network lines
for f = 1:length(network_names)
    %delete previous version
    if isfield(handles,strcat('panel_',network_names{f}))
        to_del = get(handles.(strcat('panel_',network_names{f})),'Children');
        for i = 1:1:length(to_del)
            if iscell(to_del)
                delete(to_del{i});
            else
                delete(to_del(i));
            end
        end
    end
    handles.(strcat('panel_',network_names{f})) = uipanel('Title',network_names{f},...
        'Units','characters','Tag',strcat('panel_',network_names{f}),'FontSize',12,...
        'BackgroundColor',[.94 .94 .94],'Position',[0,0,102,48],'Parent',handles.uipanelMain1,'Visible','off');
    handles.(strcat('axes_',network_names{f})) = axes('Units','characters',...
        'Position', [10,4,80,40],'Tag', strcat('axes_',network_names{f}),...
        'Parent',handles.(strcat('panel_',network_names{f})),'NextPlot','add',...
        'Xlim',[10,90],'Ylim',[4,44],'TickLength',[0,0],'XTick',[],'YTick',[],'Visible','on');
    if ~strcmp(network_names{f},'Hydro')
        loc = l_norm;
        for i = 1:1:nodes%add popup  menu to each network panel
            n_name_fix = strrep(strrep(strrep(node_names{i},' ','_'),',',''),'-','_');
            if ~isempty(network(i).(network_names{f}))
                n_e = length(network(i).Equipment);
                list={};
                list(1) = node_names(i);
                for j = 1:1:n_e
                    equip = network(i).Equipment{j};
                    r = strfind(equip,'.');
                    equip_name = equip(r+1:end);
                    k = nonzeros((1:n_g)'.*strcmp(equip_name,gen_names));
                    outs = gen_outs(gen(k(1)),nn_abrev');
                    if any(outs==f)
                        list(end+1) = {equip_name};
                    end
                end
                if length(list)>1
                    handles.(strcat('node_',n_name_fix,'_',network_names{f},'_menu')) = uicontrol('Style', 'popupmenu', 'String', list,'Units','characters','FontSize', 10,'Position', [loc(i,1)-7.5 loc(i,2)-1 15 2],'BackgroundColor','white',...
                        'Tag',strcat('node_',n_name_fix,'_',network_names{f},'_menu'),'Parent', handles.(strcat('panel_',network_names{f})),'Callback',eval(strcat('@(hObject,eventdata)MainScreen1(',quote,'node_menu',quote,',hObject,eventdata,guidata(hObject))')),'Visible','on');%popup menu with node name and all equipment
                end
            end
        end
    else
        loc = h_norm;
%         for i = 1:1:nodes
%             if ~isempty(network(i).Hydro)
%                 n_name_fix = strrep(strrep(strrep(node_names{i},' ','_'),',',''),'-','_');
%                 equip = network(i).Equipment{1};
%                 r = strfind(equip,'.');
%                 equip_name = equip(r+1:end);
%                 handles.(strcat('node_',n_name_fix)) = uicontrol('Style', 'pushbutton','String',equip_name,...
%                     'Units','characters','Position', [loc(i,1)-6 loc(i,2)-1 12 1.25],'BackgroundColor',[.85 .7 .1],...
%                     'Tag',strcat('node_',n_name_fix),'Parent',handles.(strcat('panel_',network_names{f})),...
%                     'Callback',eval(strcat('@(hObject,eventdata)MainScreen1(',quote,'edit_system',quote,',hObject,eventdata,guidata(hObject))')),'Visible','on');%button for equipment
%             end
%         end
    end
    %% add lines
    add_lines(network_names{f},subnet,node_names,loc,handles.(strcat('axes_',network_names{f})))
end
%create a panel for each node and add buttons for each piece of equipment
for i = 1:1:nodes
    n_name_fix = strrep(strrep(strrep(strrep(node_names{i},' ','_'),',',''),'-','_'),'.','');
    if isfield(handles,strcat('node_',n_name_fix))
        delete(handles.(strcat('node_',n_name_fix)));
    end
    handles.(strcat('node_',n_name_fix)) = uipanel('Title',node_names{i},...
        'Units','characters','Tag',strcat('node_',n_name_fix),'FontSize',12,...
        'BackgroundColor',[.94 .94 .94],'Position',[0,0,102,48],'Parent', handles.uipanelMain1,'Visible','off');%panel to contain equipment at this node
    nn_inc = [];%find subset of networks that have equipment at this node
    net_at_node = {};
    for f = 1:1:length(network_names)
        if ~isempty(network(i).(network_names{f}))
            nn_inc(end+1) = f;
            net_at_node(end+1) = nn_abrev(f);
        end
    end
    nn_col = nn_color(nn_inc,:);
    n_net = length(nn_inc);
    x_pos = 10 + 80*(1:n_net)'/(n_net+1);

    for k = 1:1:n_net%create bus bars
        if strcmp(network_names{nn_inc(k)},'Hydro')
            equip = network(i).Equipment{1};
            r = strfind(equip,'.');
            equip_name = equip(r+1:end);
            handles.(strcat('node_',n_name_fix,'Hydro')) = uicontrol('Style', 'pushbutton','String',equip_name,...
                'Units','characters','Position', [h_norm(i,1)-6 h_norm(i,2)-1 12 1.25],'BackgroundColor',[.85 .7 .1],...
                'Tag',strcat('node_',n_name_fix,'Hydro'),'Parent',handles.(strcat('panel_',network_names{f})),...
                'Callback',eval(strcat('@(hObject,eventdata)MainScreen1(',quote,'edit_system',quote,',hObject,eventdata,guidata(hObject))')),'Visible','on');%button for equipment
        else
            handles.(strcat('node_',n_name_fix,'_',network_names{nn_inc(k)})) = uicontrol('Style', 'pushbutton',...
                'Units','characters','Position', [x_pos(k)-1.5 4 3 40],'BackgroundColor',nn_col(k,:),...
                'Tag',strcat('node_',n_name_fix,'_',network_names{nn_inc(k)}),'Parent', handles.(strcat('node_',n_name_fix)),'Visible','on');%bus bar for energy balance
            handles.(strcat('node_',n_name_fix,'_',network_names{nn_inc(k)},'_text')) = uicontrol('Style','text',...
                'Units','characters','Position', [x_pos(k)-7.5 44 15 1.5],'BackgroundColor',[.94 .94 .94],...
                'String',network_names{nn_inc(k)},'HorizontalAlignment','center','FontSize',10,...
                'Tag',strcat('node_',n_name_fix,'_',network_names{nn_inc(k)},'_text'),'Parent', handles.(strcat('node_',n_name_fix)),'Visible','on');%bus bar for energy balance
            handles.(strcat('node_',n_name_fix,'_',network_names{nn_inc(k)},'_load')) = uicontrol('Style', 'pushbutton','String','Demand',...
                'Units','characters','Position', [x_pos(k)-6 2 12 2],'BackgroundColor',[0,0,0],'ForegroundColor',[1,1,1],...
                'Tag',strcat('node_',n_name_fix,'_',network_names{nn_inc(k)},'_load'),...
                'Callback',eval(strcat('@(hObject,eventdata)MainScreen1(',quote,'load_demands',quote,',hObject,eventdata,guidata(hObject))')),...
                'Parent', handles.(strcat('node_',n_name_fix)),'Visible','on');%bus bar for energy balance
        end
    end
    %add equipment
    n_e = length(network(i).Equipment);
    b_pos(i).y_pos = zeros(n_e,1);
    b_pos(i).x_pos = x_pos;
    b_pos(i).net = net_at_node;
    b_pos(i).nn_col = nn_col;
    for j = 1:1:n_e
        b_pos(i).y_pos(j) = 42-2*j;
        equip = network(i).Equipment{j};
        r = strfind(equip,'.');
        equip_name = equip(r+1:end);
        e_name2 = strrep(strrep(strrep(equip_name,' ','_'),',',''),'-','_');
        e_name2 = strrep(strrep(strrep(e_name2,'(',''),')',''),'.','');
        k = nonzeros((1:n_g)'.*strcmp(equip_name,gen_names));
        outs = gen_outs(gen(k(1)),net_at_node);
        if ~isempty(outs)
            handles.(strcat('node_',n_name_fix,'_',e_name2)) = uicontrol('Style', 'pushbutton','String',equip_name,...
                'Units','characters','Position', [x_pos(outs(1))-7.5 b_pos(i).y_pos(j) 15 2],'BackgroundColor',nn_col(outs(1),:),...
                'Tag',strcat('node_',n_name_fix,'_',e_name2),'Parent', handles.(strcat('node_',n_name_fix)),...
                'Callback',eval(strcat('@(hObject,eventdata)MainScreen1(',quote,'edit_system',quote,',hObject,eventdata,guidata(hObject))')),'Visible','on');%button for equipment
        end
        for se = 2:1:length(outs)%add 'buttons that are lines to other bus bars
            if outs(se)<outs(1)
                pos = [x_pos(outs(se)) b_pos(i).y_pos(j)+.85 (x_pos(outs(1))-7.5-x_pos(outs(se))) .3];
            else
                pos = [x_pos(outs(1))+7.5 b_pos(i).y_pos(j)+.85 (x_pos(outs(se))-(x_pos(outs(1))+7.5)) .3];
            end
            handles.(strcat('node_',n_name_fix,'_',e_name2,num2str(se))) = uicontrol('Style', 'pushbutton',...
                'Units','characters','Position', pos,'BackgroundColor',nn_col(outs(se),:),...
                'Tag',strcat('node_',n_name_fix,'_',e_name2,num2str(se)),'Parent', handles.(strcat('node_',n_name_fix)),'Visible','on');%line connecting to other bus bar
        end
    end
    %button to return to network view
    handles.(strcat('node_',n_name_fix,'_return')) = uicontrol('Style', 'pushbutton', 'String', 'return to network view',...
        'Units','characters','FontSize', 10,'Position', [70 .5 30 2],'BackgroundColor',[1 .4 .4],...
        'Tag',strcat('node_',n_name_fix,'_return'),'Parent', handles.(strcat('node_',n_name_fix)),...
        'Callback',eval(strcat('@(hObject,eventdata)MainScreen1(',quote,'node_view',quote,',hObject,eventdata,guidata(hObject))')),'Visible','on');%button to return to network view   
end

%add building to visualization
n_b = length(buildings);
for i = 1:1:n_b
    n_num = nonzeros((1:nodes)'.*strcmp(buildings(i).Location,node_names));
    n_name_fix = strrep(strrep(strrep(node_names{n_num},' ','_'),',',''),'-','_');
    b_name_fix = strrep(strrep(strrep(buildings(i).Name,' ','_'),',',''),'-','_');
    b_pos(n_num).y_pos(end+1) = b_pos(n_num).y_pos(end)-2;
    x_pos = b_pos(n_num).x_pos;
    outs = gen_outs(buildings(i),b_pos(n_num).net);
    handles.(strcat('node_',n_name_fix,'_',b_name_fix)) = uicontrol('Style', 'pushbutton','String',buildings(i).Name,...
                'Units','characters','Position', [b_pos(n_num).x_pos(outs(1))-12.5 b_pos(n_num).y_pos(end) 25 2],'BackgroundColor',[1,.4,.4],...%b_pos(n_num).nn_col(outs(1),:),...
                'Tag',strcat('node_',n_name_fix,'_',b_name_fix),'Parent', handles.(strcat('node_',n_name_fix)),...
                'Callback',eval(strcat('@(hObject,eventdata)MainScreen1(',quote,'edit_building',quote,',hObject,eventdata,guidata(hObject))')),'Visible','on');%button for equipment
    for se = 2:1:length(outs)%add 'buttons that are lines to other bus bars
        if outs(se)<outs(1)
            pos = [x_pos(outs(se)) b_pos(n_num).y_pos(end)+.4*se (x_pos(outs(1))-12.5-x_pos(outs(se))) .3];
        else
            pos = [x_pos(outs(1))+12.5 b_pos(n_num).y_pos(end)+.4*se (x_pos(outs(se))-(x_pos(outs(1))+12.5)) .3];
        end
        handles.(strcat('node_',n_name_fix,'_',b_name_fix,num2str(se))) = uicontrol('Style', 'pushbutton',...
            'Units','characters','Position', pos,'BackgroundColor',nn_col(outs(se),:),...
            'Tag',strcat('node_',n_name_fix,'_',b_name_fix,num2str(se)),'Parent', handles.(strcat('node_',n_name_fix)),'Visible','on');%line connecting to other bus bar
    end
end

if nodes == 1
    set(handles.(strcat('node_',n_name_fix)),'Visible','on');
    set(handles.(strcat('node_',n_name_fix,'_return')),'Visible','off');
    set(handles.network_select,'Visible','off');
else
    uistack(handles.network_select,'top')
    set(handles.network_select,'Visible','on');
    set(handles.(strcat('panel_',network_names{get(handles.network_select,'Value')})),'Visible','on');
end
end%ends function network_representation

function add_lines(net,subnet,node_names,l_norm,h)
%add new lines or remove lines from any previous network plotted
l_names = subnet.(net).lineNames;
nn = length(node_names);
l_size = 3*ones(length(l_names),1);
if isfield(subnet.(net),'lineEff')
    l_size = log(1+2*subnet.(net).lineEff);
end
for i = 1:1:length(l_names)
    %deconstruct name to find nodes that are connected
    r = strfind(l_names{i},'_');
    n1_name = l_names{i}(1:r(1)-1);
    n2_name = l_names{i}(r(2)+1:end);
    n1 = nonzeros((1:nn)'.*strcmp(n1_name,node_names));
    n2 = nonzeros((1:nn)'.*strcmp(n2_name,node_names));
    plot(h,[l_norm(n1,1);l_norm(n2,1);],[l_norm(n1,2);l_norm(n2,2);],'linewidth',l_size(i));
end
ag_nodes = length(subnet.(net).nodes);
for n = 1:1:ag_nodes
    c_nodes = subnet.(net).nodes{n};
    if length(c_nodes)>1 %interconnected perfect nodes
        n1 = nonzeros((1:nn)'.*strcmp(c_nodes{1},node_names));
%         for i = 1:1:(length(c_nodes)-1)
%             n1 = nonzeros((1:nn)'.*strcmp(c_nodes{i},node_names));
            for j = 2:1:length(c_nodes)
                n2 = nonzeros((1:nn)'.*strcmp(c_nodes{j},node_names));
                plot(h,[l_norm(n1,1);l_norm(n2,1);],[l_norm(n1,2);l_norm(n2,2);],'k','linewidth',3);
            end
%         end
    end
end
end%Ends function add_remove_lines


function out = gen_outs(gen,abbrev)
out = [];
n_net = length(abbrev);
switch gen.Type
    case 'Electric Generator'
        if isfield(gen.Output,'DirectCurrent') && any(strcmp('DC',abbrev))
            out(end+1) = nonzeros((1:n_net)'.*strcmp('DC',abbrev'));
        end
        if isfield(gen.Output,'Electricity') && any(strcmp('E',abbrev))
            out(end+1) = nonzeros((1:n_net)'.*strcmp('E',abbrev'));
        end
    case 'CHP Generator'
        if isfield(gen.Output,'DirectCurrent') && any(strcmp('DC',abbrev))
            out(end+1) = nonzeros((1:n_net)'.*strcmp('DC',abbrev'));
        end
        if isfield(gen.Output,'Electricity') && any(strcmp('E',abbrev))
            out(end+1) = nonzeros((1:n_net)'.*strcmp('E',abbrev'));
        end
        if any(strcmp('H',abbrev))
            out(end+1) = nonzeros((1:n_net)'.*strcmp('H',abbrev'));
        end
    case 'AC_DC'
            out(end+1) = nonzeros((1:n_net)'.*strcmp('E',abbrev'));
            out(end+1) = nonzeros((1:n_net)'.*strcmp('DC',abbrev'));
    case 'Utility'
        if strcmp(gen.Source,'Electricity') && any(strcmp('E',abbrev))
            out(end+1) = nonzeros((1:n_net)'.*strcmp('E',abbrev'));
        end
        if strcmp(gen.Source,'DistrictHeat') && any(strcmp('H',abbrev))
            out(end+1) = nonzeros((1:n_net)'.*strcmp('H',abbrev'));
        end
        if strcmp(gen.Source,'DistrictCool') && any(strcmp('C',abbrev))
            out(end+1) = nonzeros((1:n_net)'.*strcmp('C',abbrev'));
        end
    case 'Heater'
        if any(strcmp('H',abbrev))
            out(end+1) = nonzeros((1:n_net)'.*strcmp('H',abbrev'));
        end
    case 'Thermal Storage'
        if strcmp(gen.Source,'Heat') && any(strcmp('H',abbrev))
            out(end+1) = nonzeros((1:n_net)'.*strcmp('H',abbrev'));
        end
        if strcmp(gen.Source,'Cooling') && any(strcmp('C',abbrev))
            out(end+1) = nonzeros((1:n_net)'.*strcmp('C',abbrev'));
        end
    case 'Hydro Storage'
        if any(strcmp('W',abbrev))
            out(end+1) = nonzeros((1:n_net)'.*strcmp('W',abbrev'));
        end
        if any(strcmp('E',abbrev))
            out(end+1) = nonzeros((1:n_net)'.*strcmp('E',abbrev'));
        end
    case {'Electric Storage';'Solar'}
        if any(strcmp('DC',abbrev))
            out(end+1) = nonzeros((1:n_net)'.*strcmp('DC',abbrev'));
        elseif any(strcmp('E',abbrev))
            out(end+1) = nonzeros((1:n_net)'.*strcmp('E',abbrev'));
        end
    case 'Chiller'
        if any(strcmp('C',abbrev))
            out(end+1) = nonzeros((1:n_net)'.*strcmp('C',abbrev'));
        end
        if strcmp(gen.Source,'Heat') && any(strcmp('H',abbrev))
            out(end+1) = nonzeros((1:n_net)'.*strcmp('H',abbrev'));
        end
        if strcmp(gen.Source,'Electricity') && any(strcmp('E',abbrev))
            out(end+1) = nonzeros((1:n_net)'.*strcmp('E',abbrev'));
        end
    case 'Cooling Tower'
        if any(strcmp('CW',abbrev))
            out(end+1) = nonzeros((1:n_net)'.*strcmp('CW',abbrev'));
        end
        if any(strcmp('C',abbrev))
            out(end+1) = nonzeros((1:n_net)'.*strcmp('C',abbrev'));
        end
    case 'SingleZone'
        if any(strcmp('E',abbrev))
            out(end+1) = nonzeros((1:n_net)'.*strcmp('E',abbrev'));
        end
        if any(strcmp('H',abbrev))
            out(end+1) = nonzeros((1:n_net)'.*strcmp('H',abbrev'));
        end
        if any(strcmp('C',abbrev))
            out(end+1) = nonzeros((1:n_net)'.*strcmp('C',abbrev'));
        end
end

end%ends function gen_outs
