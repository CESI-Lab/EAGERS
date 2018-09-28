function [gen,building,fluid_loop] = manual_ic(gen,building,fluid_loop,subnet)
n_g = length(gen);
n_b = length(building);
n_fl = length(fluid_loop);
list = {};
list2 = {};
sizes = {};
index =[];
ic = zeros(1,n_g);
include = {'CHP Generator', 'Electric Generator','Chiller','Heater','Hydro Storage'};
for i = 1:1:length(gen)
    if strcmp(gen(i).Type,'Hydro Storage')
        list(end+1) = {strcat(gen(i).Name,' --- Electric Output')};
        list2(end+1) = {strcat(gen(i).Name,' --- % of max capacity')};
        sizes(end+1) = {'50'};
        index(end+1) = i;
    elseif isfield(gen(i).OpMatA,'Stor')
        list(end+1) = {strcat(gen(i).Name,' --- % of max capacity')};
        sizes(end+1) = {'50'};
        index(end+1) = i;
    elseif ismember(cellstr(gen(i).Type),include)
        list(end+1) = {strcat(gen(i).Name,' --- If greater than 0, IC must be between lower bound(',num2str(lower_bound(i)),') and upper bound(',num2str(gen(i).Size),').')};
        sizes(end+1) = {'0'};
        index(end+1) = i;
    elseif ~isempty(strcmp(gen(i).Source,'Renewable'))
        %renewable
    elseif ~isempty(strcmp(gen(i).Type,'Utility'))
        %utility
     end
end
list = [list,list2];
sizes = [sizes,50*ones(1,length(list2))];
input = str2double(inputdlg(list,'Specify Initial Condition (kW) or State of Charge (%)',1,sizes));
ic(index) = input(1:length(index));
for i=1:1:n_g
    if isfield(gen(i).OpMatA,'Stor') && ~strcmp(gen(i).Type,'Hydro Storage')
        ic(i) = ic(i)/100*gen(i).OpMatA.Stor.UsableSize; % IC = halfway charged energy storage
    end
end
%% specify initial river flow and spillway flows: IC (nL)
network_names = fieldnames(subnet);
if any(strcmp('Hydro',network_names))
    for n = 1:1:length(subnet.Hydro.nodes) 
        i = subnet.Hydro.Equipment{n};
        if strcmp(gen(i).Type,'Hydro Storage')
            gen(i).CurrentState(2) =  input(length(index)+n);
        end
    end
end
for i = 1:1:n_g
    lower_bound = 0;
    if isempty(strfind(gen(i).Type,'Storage')) && isempty(strfind(gen(i).Type,'Utility'))
        states = gen(i).OpMatB.states;
        for j = 1:1:length(states)
            lower_bound = lower_bound + gen(i).OpMatB.(states{j}).lb;
        end
    end
    gen(i).CurrentState(1) = ic(i);
    gen(i).Status = ic(i)>=lower_bound;
end
for i = 1:1:n_b
    building(i).Tzone = 20;
    building(i).Twall = 20;
    building(i).Timestamp = 0;
end
for i = 1:1:n_fl
    fluid_loop(i).fluid_temperature = 29.44; %85F
end
end%Ends function manual_ic