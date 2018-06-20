function [new_gen,new_costs,network] = design_resize(gen,equip_costs,network,size,index_size)
% for each generator specified by the input indices, update Plant with its
% newly scaled version
for j = 1:1:length(index_size)
    i = index_size(j);
    if size(j)>0
        scale = size(j) / gen(i).Size;
        gen(i) = update_component_spec(gen(i),'size', size(j));   
        equip_costs(i).Cost = equip_costs(i).Cost * scale;
        equip_costs(i).OandM = equip_costs(i).OandM * scale;
    end
end
%eliminate generators that have been set to zero size.
n_g = length(gen);
new_gen = [];
gen_names = cell(n_g,1);
for i = 1:1:n_g
    gen_names(i) = {strcat(gen(i).Type,'.',gen(i).Name)};
    if ~any(index_size==i) || size(index_size==i)>0
        if isempty(new_gen)
            new_gen = gen(i);
            new_costs = equip_costs(i);
        else
            new_gen(end+1) = gen(i);
            new_costs(end+1) = equip_costs(i);
        end
    end
end        
rem_gen = index_size(size==0);
for i = 1:1:length(rem_gen)
    rem_name = gen_names(rem_gen(i));
    for node = 1:1:length(network)
        if any(strcmp(rem_name,network(node).Equipment))
            network(node).Equipment = network(node).Equipment(~strcmp(rem_name,network(node).Equipment));
        end
    end
end
end%ends function design_resize