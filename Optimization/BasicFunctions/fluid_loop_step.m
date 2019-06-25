function [temperature,intermediate_temps] = fluid_loop_step(gen,equip,fluid_loop,temperature,best_dispatch,dt)
capacitance = fluid_loop.fluid_capacity*fluid_loop.fluid_capacitance; %Water capacity in kg and thermal capacitance in kJ/kg*K to get kJ/K
intermediate_temps = zeros(length(dt),1);
for t = 1:1:length(dt)
    imbalance = 0;
    for k = 1:1:length(equip)
        j = equip(k);
        if strcmp(gen(j).Type,'Chiller')
            imbalance = imbalance + best_dispatch(t,j) + chill_input(gen(j),best_dispatch(t,j));
        elseif strcmp(gen(j).Type,'Cooling Tower')
            imbalance = imbalance - best_dispatch(t,j);
        end
    end
    temperature = temperature + dt(t)/capacitance*imbalance;
    intermediate_temps(t) = temperature;
end
end%Ends function fluid_loop_step