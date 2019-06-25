function [flow,T_return,T_target,T_hot,T_cold] = water_equipment(building,profile,T_mains,T_tank,j)
we = building.water_use;
z = f_index(we.zone(j),building.zones.name);
if isempty(we.cold_supply_temperature_schedule{j})
    T_cold = T_mains;
else
    T_cold = profile.(we.cold_supply_temperature_schedule{j});
end
if isempty(T_tank)
    T_hot =  profile.(we.hot_supply_temperature_schedule{j});
else
    T_hot = T_tank;
end 
if isempty(we.target_temperature_schedule{j})
    T_target = T_mains;
else
    T_target = min(T_hot,profile.(we.target_temperature_schedule{j}));
end
den_water = 997;%kg/m^3
m_mix = we.peak_flow(j)*profile.(we.flow_schedule{j})*den_water*building.zones.multiplier(z); 
m_hot = m_mix.*(T_target- T_cold)./(T_hot - T_cold);
m_cold = m_mix-m_hot;
flow = m_hot;
T_return = T_cold;
% e_num = f_index(we.name{j},building.plant_demand_equip.name);
% loop_connect = building.plant_demand_equip.loop(e_num);
% n_s = max(1,length(profile.(we.latent_frac_schedule{j})));
% flow = zeros(n_s,length(building.plant_loop.name));  
% if any(strcmpi(building.plant_loop.type(loop_connect),'Heating'))
%     l = loop_connect(strcmpi(building.plant_loop.type(loop_connect),'Heating'));
%     flow(:,l) = m_hot*building.zones.multiplier(z);
% end
% if any(strcmpi(building.plant_loop.type(loop_connect),'Cooling'))
%     l = loop_connect(strcmpi(building.plant_loop.type(loop_connect),'Cooling'));
%     flow(:,l) = m_cold;
% end
% % drain_flow = m_target - m_dot_evap;
% % drain_temp = (m_target.*Cp_water.*T_target - Q_sensible - Q_latent)./(drain_flow*Cp_water);
end%Ends function water_equipment