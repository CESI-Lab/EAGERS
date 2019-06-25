function [fan_avail,override] = availability_managers(building,T_zone,T_set,dt,fan_avail)
%% Handles night cycle for HVAC system
override = false;
for i = 1:1:length(building.hvac.managers.name)
   switch building.hvac.managers.type{i}
        case 'NightCycle'
            l = building.hvac.managers.loop(i);
            z_i = f_index(1,building.hvac.zone_2_loop(l,:)');%zones on this loop
            T_max = T_set.cool' + building.hvac.managers.thermostat_tolerance(i);
            T_min = T_set.heat' - building.hvac.managers.thermostat_tolerance(i);
            switch building.hvac.managers.runtime_control_type{i}
                case 'FixedRunTime'
                    if any((T_zone(z_i)>T_max(z_i) | T_zone(z_i)<T_min(z_i)) & ~T_set.no_hvac(z_i)') %% outside bounds, need to override
                        on_loop = f_index(l,building.hvac.components.loop);
                        for j = 1:1:length(on_loop)
                            switch building.hvac.components.type{on_loop(j)}
                                case {'Fan:ConstantVolume';'Fan:VariableVolume';'Fan:OnOff'}  
                                    f = f_index(building.hvac.components.name{on_loop(j)},building.fans.name);
                                    if fan_avail(f) == 0
                                        fan_avail(f) = min(1,building.hvac.managers.runtime(i)/dt);
                                        override = true;%turn fan on when system is off
                                    end
                                case 'AirLoopHVAC:UnitaryHeatCool'
                                    u = f_index(building.hvac.components.name{on_loop(j)},building.unitary_heat_cool.name);
                                    f = f_index(building.unitary_heat_cool.fan_name{u},building.fans.name);
                                    if fan_avail(f) == 0
                                        fan_avail(f) = min(1,building.hvac.managers.runtime(i)/dt);
                                        override = true;%turn fan on when system is off
                                    end
                            end
                        end
                    end
                otherwise
                    disp('need runtime control other than fixedruntime for night manager')
            end

        case 'Scheduled'
            %not sure wat to change, if it is already scheduled, it is
            %using the profile already (except maybe air loop flows)
%             disp('need to code scheduled availability manager')
        otherwise
            disp('need to code other types of availability managers')
    end
end
end%Ends function availability_managers