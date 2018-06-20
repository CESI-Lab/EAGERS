function network = create_new_node(network,type_name)
nodes = length(network);
node_names = cell(nodes,1);
network_names = fieldnames(network);
network_names = network_names(~strcmp('name',network_names));
network_names = network_names(~strcmp('Equipment',network_names));
network_names = network_names(~strcmp('Location',network_names));
for i = 1:1:nodes
    node_names(i) = {network(i).name};
    location = zeros(nodes,3); %lat, long, time zone
    if isfield(network,'Location') && ~isempty(network(i).Location)
        location(i,1) = network(i).Location.Longitude;
        location(i,2) = network(i).Location.Latitude;
        location(i,3) = network(i).Location.TimeZone;
    else
        location(i,1:3) = nan;
    end
end
if all(isnan(location))
    default_ans = {strcat('node',num2str(nodes+1));num2str(-117.1817);num2str(46.7298);num2str(-8);};
else
    default_ans = {strcat('node',num2str(nodes+1));num2str(mean(location(~isnan(location(:,2)),2)));num2str(mean(location(~isnan(location(:,1)),1)));num2str(mean(location(~isnan(location(:,3)),3)));};
end
A = (inputdlg({'Node Name';'Location (Latitude)';'Location (Longitude)';'Location (time zone)';},'Add new node to network',1,default_ans));

network(nodes+1).name = A{1};
network(nodes+1).Location.Latitude = str2double(A{2});
network(nodes+1).Location.Longitude = str2double(A{3});
network(nodes+1).Location.TimeZone = str2double(A{4});
network(nodes+1).Equipment = {type_name};
for f = 1:1:length(network_names)
    list_nodes = {};
    net = network_names{f};
    for i = 1:1:nodes
        if ~isempty(network(i).(net))
            list_nodes(end+1) = {network(i).name};
        end
    end
    if strcmp(net,'Hydro')
        [sel,OK] = listdlg('PromptString','Downriver node', 'SelectionMode','single','ListString',list_nodes);
        if OK
            network(nodes+1).(net).connections = list_nodes(sel(1));
        end
        [sel,OK] = listdlg('PromptString','Upriver node', 'SelectionMode','multiple','ListString',list_nodes);
        if OK
            for j = 1:1:length(sel)
                k = nonzeros((1:nodes)'.*strcmp(list_nodes{sel(j)},node_names));
                network(k).(net).connections = list_nodes(sel(j));
            end
        end
        A = inputdlg({'Time to reach sea (hrs)';'Instream flow requirement (1000 cfs)';},'River segment properties',1,{'1';'0';});
        network(nodes+1).(net).Time2Sea = str2double(A{1});
        network(nodes+1).(net).InstreamFlow = str2double(A{2});
    else
        [sel,OK] = listdlg('PromptString',strcat(net,' _ connections'), 'SelectionMode','multiple','ListString',list_nodes);
        if OK
            network(nodes+1).(net).connections = cell(1,length(sel));
            network(nodes+1).(net).Trans_Eff = ones(1,length(sel));
            network(nodes+1).(net).Trans_Limit = ones(1,length(sel));
            for j = 1:1:length(sel)
                k = nonzeros((1:nodes)'.*strcmp(list_nodes{sel(j)},node_names));
                A = inputdlg({'Forward Transmission Efficiency(%)';'Forward Transmission Limit (kW)';'Reverse Transmission Efficiency(%)';'Reverse Transmission Limit (kW)';},strcat('Line properties from ',network(nodes+1).name,'_to_',network(k).name),1,{'100';'Inf';'100';'Inf';});
                network(nodes+1).(net).connections(j) = list_nodes(sel(j));
                if isempty(network(k).(net).connections)
                    network(k).(net).connections = {network(nodes+1).name};
                else
                    network(k).(net).connections(end+1) = {network(nodes+1).name};
                end
                network(nodes+1).(net).Trans_Eff(j) = str2double(A{1})/100;
                network(k).(net).Trans_Eff(end+1) = str2double(A{3})/100;
                network(nodes+1).(net).Trans_Limit(j) = str2double(A{2});
                network(k).(net).Trans_Limit(end+1) = str2double(A{4});
            end
        end
    end 
end