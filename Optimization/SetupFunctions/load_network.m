function [subnet,gen,buildings] = load_network(network,gen,buildings)
%identify generators and lines and their position in the network
%organize into sub-networks
% Group nodes between which transmission losses don't occur
% only create lines where transmission losses do occur
%% Network names, their corresponding abbreviation, what they represent
% Electrical --- E  --- standard 480V AC electrical network 
% DistrictHeat --- H --- standard 80C supply heating bus
% DistrictCool --- C --- standard 4C supply cooling bus
% Hydro        --- W --- River network with reservoirs and dams
% DirectCurrent --- DC --- 48V DC electrical network
% CoolingWater --- CW --- Water circulated between chillers and cooling towers
% Transmission1 --- E1 --- 230kV electric transmission (E2, E3, etc can be additional voltage levels
% Hydrogen     --- Hy --- Gaseous hydrogen stream
% LiqHydrogen  --- LH2 --- Liquid hydrogen
% Heating2     --- H2 --- Heat a different temperature than DistrictHeat (H3, H4... as needed)
%%-----%%%
nodes = length(network);
node_names = cell(nodes,1);
nn_list = {'Electrical';'DistrictHeat';'DistrictCool';'Hydro';'DirectCurrent';'CoolingWater';'Transmission1';'Hydrogen';'LiqHydrogen';'Heating2';};
nn_abrev = {'E';'H';'C';'W';'DC';'CW';'E1';'Hy';'LH2';'H2';};
network_names = fieldnames(network);
network_names = network_names(~strcmp('name',network_names));
network_names = network_names(~strcmp('Equipment',network_names));
network_names = network_names(~strcmp('Location',network_names));
network_names = network_names(~strcmp('Hydro',network_names));
for i = 1:1:nodes
    node_names(i) = {network(i).name};
    if isfield(network,'Location') && ~isempty(network(i).Location)
        location(i).Longitude = network(i).Location.Longitude;
        location(i).Latitude = network(i).Location.Latitude;
        location(i).TimeZone = network(i).Location.TimeZone;
    else
        location(i) = {[]};
    end
end
line_number = 0; %cumulative line segment #
for net = 1:1:length(network_names)
    subnet.(network_names{net}).nodes = {};
    n_index = nonzeros((1:length(nn_list))'.*strcmp(network_names{net},nn_list));
    subnet.(network_names{net}).abbreviation = nn_abrev{n_index};
    subnet.(network_names{net}).lineNames = {};
    subnet.(network_names{net}).lineNumber = [];
    subnet.(network_names{net}).lineLimit = [];
    subnet.(network_names{net}).lineEff = [];
    n = 0;
    for i = 1:1:nodes
        if ~isempty(network(i).(network_names{net}))
            %first check and see if this node is already part of a subNet node
            %nodes with perfect transmission are agregated into the first node in the nameList that they have perfect 2-way connection with
            [I,aNodes,connect] = agregated_node(network,node_names{i},network_names{net});
            if I == i
                n = n+1;%add a new subnet node
                subnet.(network_names{net}).nodes(n) = {aNodes};
                subnet.(network_names{net}).Location(n) = location(i);
                L = [];
                for j = 1:1:length(aNodes)
                    I = find(strcmp(aNodes{j},node_names),1,'first');
                    if isfield(network(I).(network_names{net}),'Load') && ~isempty(network(I).(network_names{net}).Load)%%note if there is a demand at this node
                        L(1,end+1) = network(I).(network_names{net}).Load;
                    end
                end
                subnet.(network_names{net}).Load(n) = {L};
                connected_nodes = {};
                for j=1:1:length(connect(:,1))
                    if ~any(strcmp(connect{j,2},aNodes))%imperfect transmission, need a line
                        [J, cNodes,~] = agregated_node(network,connect{j,2},network_names{net});
                        pconnected = node_names{J};%name of node that the connected node will be agregated into if it is perfectly connected to any others
                        connected_nodes(end+1) = {pconnected};
                        if J>i %new line connection, otherwise this was handled previously in the reverse direction
                            [eff, limit,dir] = line_prop(network,subnet.(network_names{net}).nodes{n},cNodes,network_names{net});%find forward & reverse transmission efficiency & limit
                            if strcmp(dir,'none') %no transmission (zero efficiency)
                                %do nothing
                            else
                                line_number = line_number+1;
                                if strcmp(dir,'reverse')
                                    subnet.(network_names{net}).lineNames(end+1,1) = (strcat(pconnected,'_',network_names{net},'_',node_names(i)));
                                else
                                    subnet.(network_names{net}).lineNames(end+1,1) = (strcat(node_names(i),'_',network_names{net},'_',pconnected));
                                end
                                if strcmp(dir,'dual')
                                    subnet.(network_names{net}).lineEff(end+1,1:2) = eff;
                                    subnet.(network_names{net}).lineLimit(end+1,1:2) = limit;
                                else
                                    subnet.(network_names{net}).lineEff(end+1,1) = eff;
                                    subnet.(network_names{net}).lineLimit(end+1,1) = limit;
                                end
                                subnet.(network_names{net}).lineNumber(end+1,1) = line_number;
                            end
                        end
                    end
                end
                subnet.(network_names{net}).connections{n} = connected_nodes;
            end
        end
    end
end
[subnet,gen] = create_equip_lists(network,subnet,gen,node_names);
if isfield(network,'Hydro')
    [subnet,gen] = add_river_segments(network,gen,subnet,line_number,location,node_names);
end
if ~isempty(buildings)
    [buildings,subnet] = add_buildings(buildings,subnet,gen);
end
end%Ends function load_network

function [subnet,gen] = create_equip_lists(network,subnet,gen,node_names)
%identify equipment at each subNet node (equipment can apear in multiple
%sub-nets if it produces heat and power, or uses water to produce electricity
n_g = length(gen);
gen_names = cell(n_g,1);
for i = 1:1:n_g
    gen_names(i,1) = {gen(i).Name};
end
network_names = fieldnames(subnet);
for i = 1:1:length(network_names)
    net = network_names{i};
    for n = 1:1:length(subnet.(net).nodes)
        gen_at_node = [];
        node_m = subnet.(net).nodes{n};
        for k = 1:1:length(node_m)
            net_i = nonzeros((1:length(node_names))'.*(strcmp(node_m{k},node_names)));
            equip = network(net_i).Equipment;
            for j = 1:1:length(equip)
                s = strfind(equip{j},'.');
                gen_i = nonzeros((1:length(gen_names))'.*strcmp(equip{j}(s+1:end),gen_names));
                if ~isempty(gen_i)
                    if isfield(gen(gen_i).QPform.output,subnet.(net).abbreviation)
                        gen_at_node(end+1) = gen_i;
                        gen(gen_i).QPform.(net).subnetNode = n;
                    end
                else
                    disp(strcat('generator removed from library',equip{j}))
                end
            end
        end
        subnet.(net).Equipment{n} = gen_at_node;
    end
end
end%Ends function create_equip_lists

function [trans_eff,trans_limit,dir] = line_prop(network,node1,node2,net)
%find the transmission efficiency and limit between 2 connected nodes
%if one of the nodes is perfectly connected to another node there may be
%more than one pathway connecting them, so agregate the lines
nodes = length(network);
node_names = cell(nodes,1);
for i = 1:1:nodes
    node_names(i) = {network(i).name};
end
trans_eff = zeros(1,2);
trans_limit = zeros(1,2);
node1 = unique(node1);
node2 = unique(node2);
for j = 1:1:length(node1)
    name_index = find(strcmp(node1{j},node_names),1,'first');
    for k = 1:1:length(node2)
        J = find(strcmp(node2{k},node_names),1,'first');
        %forward direction efficieny & limit
        c = find(strcmp(node2{k},network(name_index).(net).connections));
        if ~isempty(c)
            if isinf(network(name_index).(net).Trans_Limit(c)) || trans_limit(1,1)==0
                trans_eff(1,1) = network(name_index).(net).Trans_Eff(c);
                trans_limit(1,1) = network(name_index).(net).Trans_Limit(c);
            else
                weight = network(name_index).(net).Trans_Limit(c)/(trans_limit(1,1) + network(name_index).(net).Trans_Limit(c));
                trans_eff(1,1) = (1-weight)*trans_eff(1,1) + weight*network(name_index).(net).Trans_Eff(c);%weighted efficiency of 2 concurrent lines
                trans_limit(1,1) = trans_limit(1,1) + network(name_index).(net).Trans_Limit(c);
            end
        end
        %reverse direction efficiency and limit
        c = find(strcmp(node1{j},network(J).(net).connections));
        if ~isempty(c)
            if isinf(network(J).(net).Trans_Limit(c)) || trans_limit(1,2)==0
                trans_eff(1,2) = network(J).(net).Trans_Eff(c);
                trans_limit(1,2) = network(J).(net).Trans_Limit(c);
            else
                weight = network(J).(net).Trans_Limit(c)/(trans_limit(1,2) + network(J).(net).Trans_Limit(c));
                trans_eff(1,2) = (1-weight)*trans_eff(1,2) + weight*network(J).(net).Trans_Eff(c);
                trans_limit(1,2) = trans_limit(1,2) + network(J).(net).Trans_Limit(c);
            end
        end
    end
end 
if trans_eff(1,1)==0 && trans_eff(1,2)>0
    dir = 'reverse';
    trans_eff = trans_eff(1,2);
    trans_limit = trans_limit(1,2);
elseif trans_eff(1,1)>0 && trans_eff(1,2)==0
    dir = 'forward';
    trans_eff = trans_eff(1,1);
    trans_limit = trans_limit(1,1);
elseif trans_eff(1,1)==0 && trans_eff(1,2)==0
    dir = 'none';
    trans_eff = [];
    trans_limit = [];
else
    dir = 'dual';
end
end%Ends function line_prop

function [name_index,a_nodes,connect] = agregated_node(network,node,net)
%Any connected nodes with perfect bi-directional transfer are agregated into the node earliest in the list of names
%This function finds which node in the list that is
nodes = length(network);
nodeNames = cell(nodes,1);
for i = 1:1:nodes
    nodeNames(i) = {network(i).name};
end
name_index = find(strcmp(node,nodeNames),1,'first');
a_nodes = {node};
connect = cell(length(network(name_index).(net).connections),2);
connect(:,2) = network(name_index).(net).connections;
connect(:,1) = {node};
[m,~] = size(connect);
k = 0;
while k<m
    k = k+1;
    if ~any(strcmp(connect(k,2),a_nodes))%avoid looking at connections to nodes already in agregated node
        [TransEff, ~,~] = line_prop(network,connect(k,1),connect(k,2),net);
        if ~isempty(TransEff) && length(TransEff) == 2 && min(TransEff)==1 && ~strcmp(net,'Hydro')%perfect bi-directional energy transfer, hydro lines are river segments, can't agregate
            J = find(strcmp(connect{k,2},nodeNames),1,'first');
            a_nodes(end+1) = nodeNames(J);%add to list of agregated nodes
            %%add additional connections to check
            c = length(network(J).(net).connections);
            connect(end+1:end+c,1) = connect(k,2);
            connect(end-c+1:end,2) = network(J).(net).connections;
            name_index = min(J,name_index); %keep lowest number (index in list of node names
        end
        [m,~] = size(connect);
    end
end
end%Ends function agregated_node


function [subnet,gen] = add_river_segments(network,gen,subnet,line_number,location,node_names)
nodes = length(network);
subnet.Hydro.nodes = {};
subnet.Hydro.abbreviation = 'W';
subnet.Hydro.lineNames = {};
subnet.Hydro.lineNumber = [];
subnet.Hydro.lineLimit = [];
subnet.Hydro.lineMinimum = []; 
subnet.Hydro.lineTime = []; 
n_g = length(gen);
gen_names = cell(n_g,1);
for i = 1:1:n_g
    gen_names(i,1) = {gen(i).Name};
end
n = 0;
hydro_node = cell(nodes,1);
downriver = cell(nodes,1);
for i = 1:1:nodes
    if ~isempty(network(i).Hydro)
        n = n+1;%add a new subnet node
        line_number = line_number+1;   
        gen_at_node = [];
        equip = network(i).Equipment;
        for j = 1:1:length(equip)
            s = strfind(equip{j},'.');
            gen_i = nonzeros((1:length(gen_names))'.*strcmp(equip{j}(s+1:end),gen_names));
            if ~isempty(gen_i)
                if isfield(gen(gen_i).QPform.output,'W')
                    gen_at_node(end+1) = gen_i;
                    gen(gen_i).QPform.Hydro.subnetNode = n;
                    if strcmp(gen(gen_i).Type,'Hydro Storage')
                        gen(gen_i).QPform.Hydro.down_river = line_number;
                    end
                end
            end
        end        
        subnet.Hydro.nodes(n) = {{network(i).name}};
        subnet.Hydro.Location(n) = location(i);
        subnet.Hydro.connections(n) = {network(i).Hydro.connections};
        subnet.Hydro.Load(n) = {[]};%not empty if there is a water diversion for irrigation etc
        subnet.Hydro.lineNumber(n,1) = line_number;
        subnet.Hydro.lineMinimum(n,1) = 0;%update minimum instream flow later: network(i).Hydro.InstreamFlow;
        subnet.Hydro.lineLimit(n,1) = inf;     
        subnet.Hydro.Equipment{n} = gen_at_node;
        hydro_node(n) = {network(i).name}; %node names of the upriver node (origin of line segment)
        if ~isempty(network(i).Hydro.connections)
            dr_node = nonzeros((1:nodes)'.*strcmp(network(i).Hydro.connections,node_names));
            downriver(n) = {network(dr_node).name}; %node names of the downriver node (end of line segment)
            subnet.Hydro.lineNames(n,1) = {strcat(network(i).name,'_Hydro_',network(dr_node).name)};
            subnet.Hydro.lineTime(n,1) = network(i).Hydro.Time2Sea - network(dr_node).Hydro.Time2Sea; %transit time from current river to downstream river
        else
            subnet.Hydro.lineNames(n,1) = {strcat(network(i).name,'_Hydro_')};%last dam before the sea
            subnet.Hydro.lineTime(n,1) = network(i).Hydro.Time2Sea; %no downstream dam
        end
    end
end
%following indexing helps locate upstream nodes during load_matrices, and forecast_hydro
subnet.Hydro.up_river = cell(n,1);
for i = 1:1:n
    subnet.Hydro.up_river(i) = {nonzeros((1:n)'.*strcmp(hydro_node{i},downriver(1:n)))};%finds which nodes the current node is downstream of
end
end%Ends function add_river_segments

function [buildings,subnet] = add_buildings(buildings,subnet,gen)
%identify what subnet node of electrical heating and cooling the building is on, and if there are chillers/heaters connected, disengage the internal HVAC
n_b = length(buildings);
build_location = cell(n_b,1);
for i = 1:1:n_b
    build_location(i,1) = {buildings(i).Location};
end
network_names = {'Electrical';'DistrictHeat';'DistrictCool';};
for nn = 1:1:length(network_names)
    net = network_names{nn};
    for n = 1:1:length(subnet.(net).nodes)
        build = [];
        node_m = subnet.(net).nodes{n};
        for k = 1:1:length(node_m)
            for i = 1:1:n_b
                if any(strcmp(build_location{i,1},node_m(k)))
                    build(end+1) = i;
                    if strcmp(net,'Electrical')
                        buildings(i).QPform.Electrical.subnetNode = n;
                        buildings(i).QPform.Location = subnet.Electrical.Location(n);
                    elseif strcmp(net,'DistrictHeat')
                        buildings(i).QPform.DistrictHeat.subnetNode = n;
                        if ~isempty(subnet.DistrictHeat.Equipment{n})
                            equip = subnet.DistrictHeat.Equipment{n};
                            for j = 1:1:length(equip)
                                if any(strcmp(gen(equip(j)).Type,{'Heater';'CHP Generator';}))
                                    buildings(i).QPform.Heating = true;
                                end
                            end
                        end
                        if ~isempty(subnet.DistrictHeat.connections{n}) %connected to heaters at a different node
                            buildings(i).QPform.Heating = true;
                        end
                        if buildings(i).QPform.Heating
                            buildings(i).QPform.H2E = buildings(i).VariableStruct.FanPower/(1.025*(50-15))/1.225; %Flow = Heating/(Cp_Air*(Tsupply(t)-Tmix)); Power = FanPower*Flow/air density
                        end
                    elseif strcmp(net,'DistrictCool')
                        buildings(i).QPform.DistrictCool.subnetNode = n;    
                        if ~isempty(subnet.DistrictCool.Equipment{n})
                            equip = subnet.DistrictCool.Equipment{n};
                            for j = 1:1:length(equip)
                                if strcmp(gen(equip(j)).Type,'Chiller')
                                    buildings(i).QPform.Cooling = true;
                                end
                            end
                        end
                        if ~isempty(subnet.DistrictCool.connections{n}) %connected to heaters at a different node
                            buildings(i).QPform.Cooling = true;
                        end
                        if buildings(i).QPform.Cooling
                            buildings(i).QPform.C2E = buildings(i).VariableStruct.FanPower/(1.025*(25-12))/1.225; %Flow = Cooling/(Cp_Air*(Tmix-Tsupply(t))); Power = FanPower*Flow/air density
                        end
                    end
                end
            end
        end
        subnet.(net).Buildings{n} = build;
    end
end
end%Ends function add_buildings