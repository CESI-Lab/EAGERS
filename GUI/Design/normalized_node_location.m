function [l_norm,node_names] = normalized_node_location(network)
nodes = length(network);
if nodes ==1
    l_norm = [50,24];
    node_names = {network.name};
else
    location = zeros(nodes,2); %lat and long
    node_names = cell(nodes,1);
    for i = 1:1:nodes
        node_names(i) = {network(i).name};
        if isfield(network,'Location') && ~isempty(network(i).Location)
            location(i,1) = network(i).Location.Longitude;
            location(i,2) = network(i).Location.Latitude;
        else
            location(i,1:2) = nan;
        end
    end
    if (max(location(:,1))-min(location(:,1))) == 0
        l_norm = 50*ones(nodes,1);
    else
        l_norm = 10+80*(location(:,1)-min(location(:,1)))/(max(location(:,1))-min(location(:,1))); %longitude normalized to beween 10 and 90
    end
    if (max(location(:,2))-min(location(:,2)))== 0
        l_norm(:,end+1) = 24;
    else
        l_norm(:,end+1) = 4+40*(location(:,2)-min(location(:,2)))/(max(location(:,2))-min(location(:,2))); %latitude normalized to beween 4 and 44
    end
    top_row = nnz(isnan(location(:,1)));
    r = 1;
    for i = 1:1:nodes
        if isnan(location(i,1))
            l_norm(i,1) =  15 + 75*r/(top_row+1);
            l_norm(i,2) = 44;
            r = r+1;
        end
    end
end
end%Ends function normalized_node_location