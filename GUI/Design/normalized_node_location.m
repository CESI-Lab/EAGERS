function [l_norm,node_names,h_norm] = normalized_node_location(network)
nodes = length(network);
if nodes ==1
    l_norm = [50,24];
    s_norm = [50,24];
    h_norm = [50,24];
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
        l_norm(:,2) = 24;
    else
        l_norm(:,2) = 4+40*(location(:,2)-min(location(:,2)))/(max(location(:,2))-min(location(:,2))); %latitude normalized to beween 4 and 44
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
%     %% create sorted hirarchy
%     n = 8;
%     s_norm = l_norm;
%     lat_group = cell(1,n);
%     for i = 1:1:n
%         lat_group(i) = {nonzeros((1:nodes)'.*(l_norm(:,1)>100/n*(i-1) & l_norm(:,1)<=100/n*i));};
%         nl_group(i) = length(lat_group{i});
%         c_group(i) = 100/n*(i-.5);
%     end
%     [n_max,k] = max(nl_group);
%     c = 20+80*nnz(l_norm(:,1)>55)/(nnz(l_norm(:,1)>55)+nnz(l_norm(:,1)<=45));
%     center = nonzeros((1:nodes)'.*(l_norm(:,1)>45 & l_norm(:,1)<=55));
%     v_g = linspace(4,44,n_max)';%equally space vertically
%     s_norm(center,1) = v_g;
%     s_norm(center,2) = c +.25*(s_norm(center,2)-50);%bring into more vertical alignment
    
    %% create sorted hirarchy (hydropower)
    if isfield(network,'Hydro')
        down_river = {};
        h_name = {};
        h_node = [];
        to_sea = [];
        for i = 1:1:nodes
            if ~isempty(network(i).Hydro)
                h_name(end+1,1) = {network(i).name};
                h_node(end+1,1) = i;
                if isempty(network(i).Hydro.connections)
                    down_river(end+1,1) = {''};
                else
                    down_river(end+1,1) = network(i).Hydro.connections;
                end
                if isempty(network(i).Hydro.connections)
                    to_sea(end+1,1) = length(h_name);
                end
            end
        end
        j = 1;
        cc = 1;
        h_sort(1,j) = to_sea(1);
        while nnz(h_sort)<length(h_name)
            [h_sort,j] = h_branches(h_sort,h_name,down_river,to_sea,j,cc);
        end
    end
    %squeeze left to right
    keep = false(1,length(h_sort));
    for i = length(h_sort(1,:)):-1:2
        if ~any((h_sort(:,i).*h_sort(:,i-1))>0)
            h_sort(:,i-1) = h_sort(:,i-1) + h_sort(:,i);
        else
            keep(i) = true;
        end
    end
    keep(1) = true;
    h_sort = h_sort(:,keep);
    [m,n] = size(h_sort);
    h_norm = zeros(nodes,2);
    for i = 1:1:m
        for j = 1:1:n
            if h_sort(i,j)>0
                h_norm(h_node(h_sort(i,j),1),1) = 10+ 80/n*(j-.5);%10 + 80/(n-1)*(j-1);
                h_norm(h_node(h_sort(i,j),1),2) = 4+40/m*(i-.5);%4 + 40/(m-1)*(i-1);
            end
        end
    end
end
end%Ends function normalized_node_location

function [h_sort,j] = h_branches(h_sort,h_name,down_river,to_sea,j,cc)
if length(h_sort(1,:))<j
    cc = cc+1;
    h_sort(1,j) = to_sea(cc);
end
k = nonzeros((1:length(h_sort(:,j)))'.*(h_sort(:,j)>0));
k = k(end);
upr = nonzeros((1:length(h_name))'.*strcmp(h_name{h_sort(k,j)},down_river));%find upriver
if ~isempty(upr)%follow each branch until it doesn't have an upstream dam 
    h_sort(k+1,j) = upr(1);%if multiple upriver, create branches
    if length(upr)>1%if multiple upriver, create branches
        n_branches = length(upr)-1;
        m_branches = length(h_sort(1,:))-j;%branches to shift right
        if m_branches>0
            h_sort(:,j+n_branches+1:j+n_branches+m_branches) = h_sort(:,j+1:j+m_branches);
            h_sort(:,j+1:j+n_branches) = 0;
        end
        h_sort(k+1,j+1:j+n_branches) = upr(2:end)';
    end
else
    j = j+1;%move to next branch
end
end%Ends function h_branches