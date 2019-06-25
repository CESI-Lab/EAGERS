function setpoint = manager_override(building,profile,node,setpoint,air,fresh_air,fan_avail)
%% overides temperature setpoints at specific nodes, e.g. cold air supply
k = f_index(node,building.manager.node);
if ~isempty(k)
    for j = 1:1:length(k)
        i = k(j);%
        switch building.manager.type{i}
            case 'OutdoorAirReset'
                if building.manager.low_temp(i)<building.manager.high_temp(i)
                    if fresh_air.T < building.manager.low_temp(i)
                        setpoint = building.manager.setpoint_at_low_temp(i);
                    elseif fresh_air.T>building.manager.high_temp(i)
                        setpoint = building.manager.setpoint_at_high_temp(i);
                    else
                        setpoint = building.manager.setpoint_at_low_temp(i) - (fresh_air.T - building.manager.low_temp(i))/(building.manager.high_temp(i) - building.manager.low_temp(i))*(building.manager.setpoint_at_low_temp(i) - building.manager.setpoint_at_high_temp(i));
                    end
                else
                    setpoint = 0.5*(building.manager.setpoint_at_low_temp(i) + building.manager.setpoint_at_high_temp(i));
                end
            case 'Scheduled'
                setpoint = profile.(building.manager.schedule{i});
            case 'MixedAir'
                if air.m_dot>0
                    f = f_index(building.manager.fan_inlet_node{i},building.fans.inlet);
                    [~,~,Q_to_air] = fan_calc(building.fans,building.fans.name{f},air,fan_avail(f));
                    %cooling in zone calculation is m*Cp*dT, not h2-h1, so need to ensure sensible load = m*Cp*(Tsupply -Tzone)
                    spec_heat_supply = 1006 + 1860*air.w; %J/kg*K
                    setpoint = setpoint - Q_to_air./(air.m_dot.*spec_heat_supply);
                end
            case 'FollowOutdoorAirTemperature'
                switch building.manager.ref_temp_type{i}
                    case {'OutdoorAirWetBulb';'OutdoorWetBulb';}
                        setpoint = psychometric(fresh_air,'Twb') + building.manager.offset_temperature(i);
                    case {'OutdoorAirDryBulb';'OutdoorDryBulb';}
                        setpoint = fresh_air.T + building.manager.offset_temperature(i);
                end
                setpoint = max(building.manager.min_temperature(i),min(building.manager.max_temperature(i),setpoint));
            case 'SingleZone:Reheat'
                %handled elswhere
            case 'SingleZone:Humidity:Minimum'
                %need to add for hospital see page 1396
            case 'SingleZone:Humidity:Maximum'
                %need to add for hospital
            otherwise
                disp('need more HVAC manager types')
                %do nothing
        end
    end
end
end%Ends function manager_override