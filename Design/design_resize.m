function [new_gen,new_costs,network] = design_resize(gen,equip_costs,network,size,multiple,index_size)
% for each generator specified by the input indices, update Plant with its
% newly scaled version
n_g = length(gen);
rem_gen = [];
for j = 1:1:length(index_size)
    i = index_size(j);
    if size(j)==0 || multiple(j)==1 %remove generator
        rem_name = strcat(gen(i).Type,'.',gen(i).Name);
        for node = 1:1:length(network)
            if any(strcmp(rem_name,network(node).Equipment))
                network(node).Equipment = network(node).Equipment(~strcmp(rem_name,network(node).Equipment));
            end
        end
        rem_gen(end+1) = i;
    end
    if size(j)>0 && size(j)~=gen(i).Size && multiple(j)>=1%change the size of the specified generator
        scale = size(j) / gen(i).Size;
        gen(i) = update_component_spec(gen(i),'size', size(j));   
        equip_costs(i).Cost = equip_costs(i).Cost * scale;
        equip_costs(i).OandM = equip_costs(i).OandM * scale;
    end
    if size(j)>0 && multiple(j) >1 %add multiples of the current generator
        primary_name = {strcat(gen(i).Type,'.',gen(i).Name)};
        for k = 2:1:multiple(j)
            n_g = n_g+1;
            gen(n_g) = gen(i);
            gen(n_g).Name = strcat(gen(n_g).Name,'_',num2str(k));
            equip_costs(n_g) = equip_costs(i);
        end
        for node = 1:1:length(network)
            if any(strcmp(primary_name,network(node).Equipment))
                for k = 2:1:multiple(j)
                    network(node).Equipment(end+1) = strcat(primary_name,'_',num2str(k));
                end
            end
        end
    end
end
%eliminate generators that have been set to zero size or multiplier is zero
new_gen = [];
for i = 1:1:n_g
    if ~any(i == rem_gen)
        if isempty(new_gen)
            new_gen = gen(i);
            new_costs = equip_costs(i);
        else
            new_gen(end+1) = gen(i);
            new_costs(end+1) = equip_costs(i);
        end
    end
end        

end%ends function design_resize