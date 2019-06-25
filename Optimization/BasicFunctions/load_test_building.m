function test_building = load_test_building(buildings,network,date,weather)
n_b = length(buildings);
n_s = length(date);
test_building.NonHVACelectric = zeros(n_s,n_b);
test_building.InternalGains = zeros(n_s,n_b);
nodes = length(network);
node_names = cell(nodes,1);
for i = 1:1:nodes
    node_names(i) = {network(i).name};
end
for i = 1:1:n_b
    net_index = nonzeros((1:nodes)'.*(strcmp(buildings(i).Location,node_names)));
    if isfield(network(net_index),'Location')
        location = network(net_index).Location;
    else
        location = {[]};
    end
    sg = solar_gain(buildings(i),date,location,weather);
    b_loads = building_loads(buildings(i),date,sg);
    test_building.InternalGains(:,i) = b_loads.InternalGains;
    test_building.NonHVACelectric(:,i) = b_loads.Equipment + b_loads.InteriorLighting + b_loads.ExteriorLighting + b_loads.OtherLoads;
    if ~isempty(b_loads.DCloads)
        test_building.DCloads(:,i) = b_loads.DCloads;
    end
end
