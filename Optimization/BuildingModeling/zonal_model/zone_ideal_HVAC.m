function [air_nodes,central_flow,direct_flow,mixed_air_loop,mixed_air_zone,mix,plenum,profile,T,w,fan_avail,window_gain] = ...
    zone_ideal_HVAC(building,profile,air_nodes,weather,mixing,infiltration,gains,T,w,T_set,w_set,occupancy,vect_2_sun,dt,des_day)
%% Iterate to find central air and direct zone air temperature and flow that achieves temperature constraints
air_density = 1.204; %kg/m^3 %standard air density (1.204 kg/m3) adjusted for the local barometric pressure (standard barometric pressure corrected for altitude, ASHRAE 1997 HOF pg. 6.1).
fresh_air.T = weather.DrybulbC;
fresh_air.w = w.air;
fresh_air.h = psychometric(fresh_air,'h');
cp_zone = (1006 + 1860*w.zone)*air_density;%Specific heat in J/m^3;
Q_adjustment = zeros(length(T.zone),1);
T.zone_est = T.zone; T.surf_est = T.surf; w.zone_est = w.zone;%estimated values at the end of the time step
T.central = [];
fresh_air_zone_req = design_zone_outdoor_air(building,occupancy);%requested fresh air of the zone
net_H2O_gain = gains.zone_latent'/2264705 + air_density*(mixing*w.zone + infiltration'.*w.air);%flows in m^3/s converted to kg/s, everything into kg/s of water 
fan_avail = ones(length(building.fans.name),1);
for i = 1:1:length(building.fans.name)
    fan_avail(i) = profile.(building.fans.schedule{i});
end
window_gain = window_solar_transmittance(building.windows,building.surfaces,vect_2_sun,weather,gains.window);
max_e = 10*building.site.temp_tol;
count = 0;
while count<8 && max_e>0.1*building.site.temp_tol
    [T.windows,T.windows_int,~,window_h] = window_temperature(building,vect_2_sun,T.windows,T.sky+273,T.zone_est,T.surf_est,cp_zone,weather,gains.window);
    dT = T.surf_est(building.ctf.sur_state1) - building.ctf.map_sur2zone'*T.zone_est;%interior surface convection = h*A*(Tsur - Tzone) 
    h_interior = natural_convection(building.surfaces.normal,dT,T_set.no_hvac,building.convection.interior);
    %% Predicts the zone energy exchange that achieves the temperature setpoints 
    %% The actual HVAC cooling/heating is different than zone cooling/heating due to recirculation & fresh air
    Q_mixing = cp_zone.*(mixing*T.zone_est - sum(mixing,2).*T.zone_est);%heat into zone via mixing from another zone
    Q_infiltration = cp_zone.*infiltration'.*(T.air_zone - T.zone_est);%heat into zone from air infiltration 
    Q_surfaces(:,1) = -building.ctf.map_sur2zone*(h_interior.*(building.ctf.map_sur2zone'*T_set.cool' - T.surf_est(building.ctf.sur_state1)).*building.surfaces.area);%estimated heat transfer to zone from walls/floor/ceiling
    Q_surfaces(:,2) = -building.ctf.map_sur2zone*(h_interior.*(building.ctf.map_sur2zone'*T_set.heat' - T.surf_est(building.ctf.sur_state1)).*building.surfaces.area);%estimated heat transfer to zone from walls/floor/ceiling
    Q_transient = (T.zone - [T_set.cool',T_set.heat']).*cp_zone.*building.zones.volume./dt;%zone heat capacitance*deltaT/dt = power (W)
    Q_windows = building.ctf.map_win2zone*(window_h.*building.windows.area.*(T.windows_int - T.zone_est(building.windows.zone_index)));% + windows.solar_gain;
    Q_Tset = -(gains.zone_sensible' + Q_windows + Q_mixing + Q_infiltration + Q_transient + Q_surfaces + Q_adjustment);%Energy provided to the zone in W to achieve T_set.cool, + is heating
    Q_Tset(T_set.no_hvac,:) = 0;
    %% determine the supply flow rate and temperature to each zone/loop that will result in the desired zone temperature (passive control within deadband)
    [T_min,T_max] = zone_limits(building,T.zone_est);
    exhaust_flow = exhaust_fan_flow(building,profile,mixing);
    flow_imbalance = (infiltration' - exhaust_flow - sum(mixing,1)' + sum(mixing,2));%volumetric flow (m^3/s)

    if isempty(des_day)
        central_flow = terminal_flow(building,Q_Tset,T_set,cp_zone,T_min,T_max,T.central,flow_imbalance,fan_avail);
        [min_loop_flow,max_loop_flow] = supply_flow_limit(building,profile,fan_avail);%enforces minimum and maximum flow rates
        loop_flow = max(min_loop_flow,min(max_loop_flow,building.hvac.zone_2_loop*(central_flow.*building.zones.multiplier)));%volumetric flow (m3/s)
        loop_2_zone = divide_loop_flow(building,loop_flow,central_flow);
        central_flow = loop_2_zone*loop_flow./building.zones.multiplier;%loop flow adjusted in air loop flow constraints were active
        direct_flow = zone_hvac_flow(building,profile,Q_Tset);
        fresh_air_loop = loop_outdoor_air(building,profile,fresh_air_zone_req,central_flow,loop_flow,flow_imbalance);
        direct_fresh_air = direct_zone_outdoor_air(building,profile,Q_Tset,flow_imbalance,central_flow);
    else
        [central_flow,direct_flow,loop_flow,direct_fresh_air,fresh_air_loop] = des_day_zone_flow(building,fresh_air_zone_req,cp_zone,Q_Tset,T_set,T_min,T_max,flow_imbalance);
    end
    [flow_exiting_zone,return_flow,mix,plenum] = plenum_return_flow(building,profile,central_flow,direct_fresh_air,mixing,infiltration');
    [direct_return_air,loop_return_air] = return_air_calc(building,gains,T.zone_est,w.zone_est,return_flow,direct_flow,direct_fresh_air);
    [T_loop,w.loop,mixed_air_loop,air_nodes] = loop_supply_setpoint(building,profile,air_nodes,fresh_air,fresh_air_loop,loop_return_air,loop_flow,central_flow,loop_2_zone,cp_zone,flow_exiting_zone,net_H2O_gain,T_set,w_set,Q_Tset,T_min,T_max,fan_avail,des_day);
    [T.central,T.direct,mixed_air_zone] = zones_supply_setpoint(building,central_flow,cp_zone,T_loop,direct_return_air,fresh_air,direct_fresh_air,direct_flow,T.zone_est,T_set,Q_Tset,T_min,T_max,fan_avail);
    
    %% Estimate of zone conditions if HVAC logic is applied
    w.central = building.hvac.zone_2_loop'*w.loop;  
    w.direct = mixed_air_zone.w;
    
    [T_z_new,T.surf_est,w.zone_est,T.windows_int] = update_model_states_implicit(building,T,gains,window_gain,dt,w,mix,plenum,infiltration',central_flow,direct_flow,weather);
    max_e = max(abs(T_z_new - T.zone_est).*(1-T_set.no_hvac'));
    %% adjustment to help converge to setpoint with temperature dynamics (not necessary?)
    Q_applied = central_flow.*cp_zone.*(T.central - [T_set.cool',T_set.heat']) + direct_flow.*cp_zone.*(T.direct - [T_set.cool',T_set.heat']);%zone energy exchange
    target_cool = abs((Q_applied(:,1) - Q_Tset(:,1))./Q_Tset(:,1))<1e-3;
    target_heat = abs((Q_applied(:,2) - Q_Tset(:,2))./Q_Tset(:,2))<1e-3;
    if any(target_cool) || any(target_heat)
        Q_ss = Q_Tset + Q_transient;
        Q_slope = (Q_ss(:,1) - Q_ss(:,2))./(T_set.cool-T_set.heat)';
        if any(target_cool)
            Q_adj_cool = (T_z_new - T_set.cool').*Q_slope;
            Q_adjustment(target_cool) = Q_adjustment(target_cool) + 0.25*Q_adj_cool(target_cool);
            max_e = max(max_e,max(abs(T_z_new(target_cool) - T_set.cool(target_cool)')));
        end
        if any(target_heat)
            Q_adj_heat = (T_z_new - T_set.heat').*Q_slope;
            Q_adjustment(target_heat) = Q_adjustment(target_heat) + 0.25*Q_adj_heat(target_heat);
            max_e = max(max_e,max(abs(T_z_new(target_heat) - T_set.heat(target_heat)')));
        end
    end
    cp_zone = (1006 + 1860*w.zone_est)*air_density;%Specific heat in J/m^3
    
    %% add night manager here, if it turns things on or off, don't update T.zone_est, and don't count up
    [fan_avail,override] = availability_managers(building,T_z_new,T_set,dt,fan_avail);
    if ~override
        T.zone_est = T_z_new;
        count = count+1;
    end
end
end%Ends function zone_ideal_HVAC

function [flow_exiting_zone,return_flow,mixing,plenum_flow] = plenum_return_flow(building,profile,central_flow,direct_fresh_air,mixing,infiltration)
%%exhaust fans
n = length(building.zones.name);
exhaust_flow = zeros(n,1);
if any(building.fans.exhaust)
    for i = 1:1:length(building.fans.name)
        if building.fans.exhaust(i) && ~isnan(building.fans.flow_rate(i))
            exhaust_flow(building.fans.exhaust_zone(i)) = building.fans.flow_rate(i)*profile.(building.fans.schedule{i});%volumetric flow (m^3/s)  %See page 1601 of the reference manual
        end
    end
end
exhaust_flow = exhaust_flow - sum(mixing,1)';%balanced exhaust flow (m^3/s)
return_flow = central_flow + direct_fresh_air + sum(mixing,2) - sum(mixing,1)' + infiltration - exhaust_flow;%volumetric flows (m3/s)
plenum_flow = zeros(length(building.zones.name));%keep track of what zone plenum flow is from
for i = 1:1:length(building.hvac.plenum.name)
    z = f_index(building.hvac.plenum.zone{i},building.zones.name);
    z2P = building.hvac.zone_2_plenum(i,:);
    z2P(z) = 0; %don't count m_sys or infiltration into a plenum zone as mixing from another zone
    plenum_flow(z,:) = z2P.*return_flow';%flow from column j to row i
end
return_flow = return_flow + sum(plenum_flow,2) - sum(plenum_flow,1)';
% mix_and_plenum = mixing+plenum_flow;%mixing flow between zones from column j to row i
flow_exiting_zone = central_flow + direct_fresh_air + sum(mixing,2) + infiltration + sum(plenum_flow,2);%flow leaves through exhaust, mixing to another zone or plenum, or return air
% loop_return_flow = central_flow + direct_fresh_air + infiltration - exhaust_flow - sum(mix_and_plenum,1)' + sum(mix_and_plenum,2);
end%Ends function plenum_return_flow

function exhaust_flow = exhaust_fan_flow(building,profile,mixing)
%%exhaust fans
n = length(building.zones.name);
exhaust_flow = zeros(n,1);
if any(building.fans.exhaust)
    for i = 1:1:length(building.fans.name)
        if building.fans.exhaust(i) && ~isnan(building.fans.flow_rate(i))
            exhaust_flow(building.fans.exhaust_zone(i)) = building.fans.flow_rate(i)*profile.(building.fans.schedule{i});%volumetric flow (m^3/s)  %See page 1601 of the reference manual
        end
    end
end
exhaust_flow = exhaust_flow - sum(mixing,1)';%balanced exhaust flow (m^3/s)
end%Ends function zone_exhaust_fans

function central_flow = terminal_flow(building,Q_Tset,T_set,cp_zone,T_min,T_max,T_central,flow_imbalance,fan_avail)
%%do VAV terminal air calculations  
%Determine flow rate to zones, based on fan flow limits, fresh air requirements and if HVAC is on
air_density = 1.204; %kg/m^3 %standard air density (1.204 kg/m3) adjusted for the local barometric pressure (standard barometric pressure corrected for altitude, ASHRAE 1997 HOF pg. 6.1).
term_flow = zeros(length(building.hvac.terminals.name),1); 
zone_flow_req = zeros(length(T_min),1);
h = Q_Tset(:,2)>0 & ~T_set.no_hvac';
zone_flow_req(h) = Q_Tset(h,2)./(cp_zone(h).*(T_max(h) - T_set.heat(h)'));
c =  Q_Tset(:,1)<0 & ~T_set.no_hvac';
zone_flow_req(c) =  Q_Tset(c,1)./(cp_zone(c).*(T_min(c) - T_set.cool(c)'));

z_list = building.hvac.terminals.zone;
min_term_flow = term_flow;
if isfield(building.hvac.terminals,'min_flow') 
    min_term_flow = building.hvac.terminals.min_flow./building.zones.multiplier(z_list);
    min_term_flow(isnan(min_term_flow)) = 0;
end
max_term_flow = building.hvac.terminals.max_flow./building.zones.multiplier(z_list);

for i = 1:1:length(building.hvac.loop.name)
    term_i = f_index(i,building.hvac.terminals.loop);
    on_loop = f_index(i,building.air_supply_equip.loop);
    k = f_index({'Fan:ConstantVolume';'Fan:VariableVolume';'Fan:OnOff';'Fan:SystemModel';},building.air_supply_equip.type(on_loop));
    if ~isempty(k)
        fan_name = building.air_supply_equip.name{on_loop(k)};
        f = f_index(fan_name,building.fans.name);
        min_term_flow(term_i) = min_term_flow(term_i)*fan_avail(f);
        max_term_flow(term_i) = max_term_flow(term_i)*fan_avail(f);
    end
    
    %Cooling mode:supply temperature is based on limiting zone, other zones operate with reduced flow. If a terminal minimum flow is hit, re-heat prevents from overcooling
    %Heating mode: minimum flow. re_heat does all necessary heating
    z_i = building.hvac.terminals.zone(term_i);
    m = building.zones.multiplier(z_i);
    switch building.hvac.terminals.type{term_i(1)}%what if terminals on same loop are of different variety? Not an issue with default buildings
        case 'Uncontrolled'
            term_flow(term_i) =  zone_flow_req(z_i);%volumetric flow m^3/s
            %check fan schedule
        case 'ConstantVolume:Reheat'
            disp('has constantvolume:reheat, dont know what to do')
%             flow_max = building.hvac.terminals.max_flow(term_i).*building.hvac.terminals.max_reheat_frac(term_i);%volumetric flow m^3/s
%             if ~any(isnan(building.hvac.terminals.max_reheat_flow_per_area(term_i)))
%                 flow_max = min(flow_max,building.hvac.terminals.max_reheat_flow_per_area(term_i).*building.zones.floor_area(z_i));%volumetric flow m^3/s
%             end
        case 'VAV:Reheat'
            if length(z_i) == 1 %serves a single zone, heating = min flow, vary temp ... cooling = supply temp, vary flow
                if c(z_i)%cooling
                    term_flow(term_i) = zone_flow_req(z_i);
                else%if not cooling, set to minimuim (bottom of page 694)
                    term_flow(term_i) = min_term_flow(term_i);
                end
            else %loop serves multiple zones, must supply at temperature that provides sufficient cooling to all zones at max flow
                %find the limiting zone, i.e. max flow at this supply temperature = perfect cooling. Other zones have re-heat
                if any(c(z_i))
                    min_cooling_T = min(max(T_min(z_i),T_set.cool(z_i)' + Q_Tset(z_i,1)./(cp_zone(z_i).*max_term_flow(term_i))));
                    if ~isempty(T_central)
                        min_cooling_T = min(min_cooling_T,min(T_central(z_i)));
                    end
                    term_flow(term_i) = Q_Tset(z_i,1)./(cp_zone(z_i).*(min_cooling_T - T_set.cool(z_i)'));%volumetric flow m^3/s
                else%passive or heating
                    flow_req = zeros(length(z_i),1);
                    for j = 1:1:length(z_i)
                        switch building.hvac.terminals.reheat_coil_type{term_i(j)}
                            case 'Coil:Heating:Water'
                                name = building.hvac.terminals.reheat_coil_name{term_i(j)};
                                c_num = f_index(name,building.coils.heating.name);
                                load_coil = min(Q_Tset(z_i(j),2)*m(j),building.coils.heating.capacity(c_num));
                                air = water_coil(building,name,air_density*term_flow(term_i(j))*m(j),load_coil,[]);
                                flow_req(j) = air.m_dot/air_density/m(j);
                            otherwise
                                flow_req(j) = Q_Tset(z_i(j),2)./(cp_zone(z_i(j)).*(T_max(z_i(j))-T_set.heat(z_i(j))));
                        end
                    end
                    term_flow(term_i) = max(flow_req,min_term_flow(term_i));%volumetric flow m^3/s
                end
            end
        otherwise
            disp('need additional air terminal types')
    end
end
term_flow =  max(term_flow, -flow_imbalance(z_list).*(~T_set.no_hvac(z_list)'));%volumetric flow m^3/s
term_flow =  min(max(term_flow,min_term_flow),max_term_flow);
central_flow = zeros(length(building.zones.name),1);
central_flow(building.hvac.terminals.zone) = term_flow;
end%Ends function terminal_flow

function direct_flow = zone_hvac_flow(building,profile,Q_Tset)
%flow for zones not on a HVAC loop (unit heaters)
direct_flow = zeros(length(building.zones.name),1);
for i = 1:1:length(building.unitary_sys.name)
    f = f_index(building.unitary_sys.fan_name{i},building.fans.name);
    z = building.unitary_sys.zone(i);
    switch building.unitary_sys.type{i}
        case 'UnitHeater'
            m = building.zones.multiplier(z);
            direct_flow(z) = building.fans.flow_rate(f)*profile.(building.fans.schedule{i})/m;
        case 'FourPipeFanCoil'
            m = building.zones.multiplier(z);
            direct_flow(z) = building.fans.flow_rate(f)*profile.(building.fans.schedule{i})/m;
        case 'PackagedTerminalAirConditioner'
            m = building.zones.multiplier(z);
            %% need to find part-load-ratio (PLR) see page 1541
            %% later infer PLR from this flow
            if Q_Tset(z,2)>0 %heating
                direct_flow(z) = building.unitary_sys.heating_air_flow(i)/m;
            elseif Q_Tset(:,1)<0 %cooling
                direct_flow(z) = building.unitary_sys.cooling_air_flow(i)/m;
            else
                direct_flow(z) = building.unitary_sys.no_load_air_flow(i)/m;
            end
        otherwise
            disp('other type of zone hvac')
    end
end
end %Ends function zone_hvac_flow

function fresh_air_loop = loop_outdoor_air(building,profile,fresh_air_zone_req,central_flow,loop_flow,flow_imbalance)
%fresh air specified by a AirLoop HVAC
nl = length(building.hvac.loop.name);
fresh_air_loop = zeros(nl,1);
for i = 1:1:length(building.hvac.outdoor_air.name)%loop through outside_air controllers to get minimum
    l = building.hvac.outdoor_air.loop(i);
    fresh_air_loop(l) = building.hvac.outdoor_air.min_flow(i)*profile.(building.hvac.outdoor_air.min_air_schedule{i});
end
fresh_air_loop = max(fresh_air_loop,-building.hvac.zone_2_loop*(flow_imbalance.*building.zones.multiplier));%mass flow that must be made up by fresh air(kg/s)
% oa_frac = fresh_air_loop./loop_flow;
% oa_frac(isinf(oa_frac) | isnan(oa_frac)) = 0;
% oa_frac_req = fresh_air_zone_req./central_flow;
% oa_frac_req(isinf(oa_frac_req) | isnan(oa_frac_req)) = 0;
% for i = 1:1:nl
%     if loop_flow(i)>0
%         z_i = f_index(1,building.hvac.zone_2_loop(i,:));
%         loop_oa_frac_req = max(oa_frac_req(z_i));
%         if loop_oa_frac_req>oa_frac(i)
%             fresh_air_loop(i) = loop_oa_frac_req*loop_flow(i);
%         end
%     end
% end
end%Ends function loop_outdoor_air

function direct_fresh_air = direct_zone_outdoor_air(building,profile,Q_Tset,flow_imbalance,central_flow)
%fresh air specified by ZoneHVAC units
direct_fresh_air = zeros(length(building.zones.name),1);
for i = 1:1:length(building.unitary_sys.name)
    z = building.unitary_sys.zone(i);
    switch building.unitary_sys.type{i}
        case 'FourPipeFanCoil'
            m = building.zones.multiplier(z);
            if ~isempty(building.unitary_sys.outdoor_air_schedule{i})
                direct_fresh_air(z) = building.unitary_sys.outdoor_air_flow(i)*profile.(building.unitary_sys.outdoor_air_schedule{i})/m;
            else
                direct_fresh_air(z) = building.unitary_sys.outdoor_air_flow(i)/m;
            end
        case 'PackagedTerminalAirConditioner'
            m = building.zones.multiplier(z);
            if Q_Tset(z,2)>0 %heating
                direct_fresh_air(z) = building.unitary_sys.heating_outdoor_air_flow(i)*profile.(building.unitary_sys.fan_schedule{i})/m;
            elseif Q_Tset(:,1)<0 %cooling
                direct_fresh_air(z) = building.unitary_sys.cooling_outdoor_air_flow(i)*profile.(building.unitary_sys.fan_schedule{i})/m;
            else
                direct_fresh_air(z) = building.unitary_sys.no_load_outdoor_air_flow(i)*profile.(building.unitary_sys.fan_schedule{i})/m;
            end
        otherwise
            %no direct outdoor airflow
    end
end
direct_fresh_air = max(direct_fresh_air,-(flow_imbalance + central_flow));
end%Ends function minimum_outdoor_air

function loop_2_zone = divide_loop_flow(building,loop_flow,central_flow)
loop_2_zone = 0*building.hvac.zone_2_loop';
for i = 1:1:length(building.hvac.loop.name)
    if loop_flow(i)>0 %split flow based on heating/cooling request flow (avoid NaN from dividing by zero)
        min_flow_loop_i = building.hvac.zone_2_loop(i,:)*(central_flow.*building.zones.multiplier);
        if min_flow_loop_i>0
            loop_2_zone(:,i) = (central_flow.*building.hvac.zone_2_loop(i,:)'.*building.zones.multiplier)/min_flow_loop_i;
        else %avoid NaN, split to zones based on max zone flow
            zone_request = 0*central_flow;
            z = f_index(1,building.hvac.zone_2_loop(i,:)');
            for j = 1:1:length(z)
                t_z = f_index(z(j),building.hvac.terminals.zone);
                zone_request(z(j)) = building.hvac.terminals.max_flow(t_z);
            end
            loop_2_zone(:,i) = zone_request/sum(zone_request);
        end
    end
end
end%Ends function divide_loop_flow

function [minimum_flow,limit_flow] = supply_flow_limit(building,profile,fan_avail)
%reduce flow if any equipment in air_loop limits flow
limit_flow = building.hvac.design.max_loop_flow;%volumetric flow (m3/s)
minimum_flow = 0*limit_flow;%volumetric flow (m3/s)
%%go through equipment and reduce if schedule is reduced
for i = 1:1:length(building.air_supply_equip.name)
    type = building.air_supply_equip.type{i};
    name = building.air_supply_equip.name{i};
    loop = building.air_supply_equip.loop(i);
    switch type
         case {'CoilSystem:Cooling:DX';'Coil:Cooling:Water'}
            k = f_index(name,building.coils.cooling.name);
            sched = building.coils.cooling.schedule{k};
            limit_flow(loop) = min(limit_flow(loop),building.coils.cooling.rated_air_flow(k)*profile.(sched));%mass flow (kg/s)
        case {'Fan:ConstantVolume';'Fan:OnOff';'Fan:VariableVolume';'Fan:SystemModel';}
            f = f_index(name,building.fans.name);
            [minimum_flow(loop),limit_flow(loop)] = fan_flow_limit(building,fan_avail(f),minimum_flow(loop),limit_flow(loop),f);               
        case 'AirLoopHVAC:UnitaryHeatCool'
            u = f_index(name,building.unitary_heat_cool.name);
            f = f_index(building.unitary_heat_cool.fan_name{u},building.fans.name);
            [minimum_flow(loop),limit_flow(loop)] = fan_flow_limit(building,fan_avail(f),minimum_flow(loop),limit_flow(loop),f);          
    end     
end
end%Ends function supply_flow_limit

function [min_flow,limit_flow] = fan_flow_limit(building,fan_avail,minimum,limit,f)
switch building.fans.type{f}
    case {'ConstantVolume';'OnOff';}
        limit_flow = building.fans.flow_rate(f)*fan_avail;%volumetric flow (m3/s)
        min_flow = building.fans.flow_rate(f)*fan_avail;%volumetric flow (m3/s)
    case 'VariableVolume'
        limit_flow = min(limit,building.fans.flow_rate(f)*fan_avail);%volumetric flow (m3/s)
        min_flow = max(minimum,building.fans.min_flow_rate(f)*fan_avail);%volumetric flow (m3/s)   
    case 'Fan:SystemModel'
        disp('need to program Fan:SystemModel')
end
end%Ends function fan_flow_limit

function [direct_return_air,loop_return_air] = return_air_calc(building,gains,T,w,return_flow,direct_flow,direct_flow_fresh)
%determine air returning to HVAC from zone, and agregate into loops
air_density = 1.204; %kg/m^3 %standard air density (1.204 kg/m3) adjusted for the local barometric pressure (standard barometric pressure corrected for altitude, ASHRAE 1997 HOF pg. 6.1).
zone_2_loop = building.hvac.zone_2_loop;
multiplier = building.zones.multiplier;

%add sensible and latent heat from refrigeration equipment
return_air.T = T;
return_air.w = w;
return_air.m_dot = return_flow*air_density;%mass flow (kg/s)
latent_H2O = gains.hvac_latent'/2264705;%convert J/s to kg/s of water
air_H2O = return_air.w.*return_air.m_dot;
return_air.m_dot = return_air.m_dot + latent_H2O;
return_air.w = (air_H2O + latent_H2O)./return_air.m_dot;
cp_zone = 1006 + 1860*return_air.w; %J/kg*K
newT = return_air.T + gains.hvac_sensible'./(cp_zone.*return_air.m_dot);
A = return_air.m_dot>0;
return_air.T(A) = newT(A);
return_air.w(~A) = w(~A);
return_air.h =  psychometric(return_air,'h');

central_return_air = return_air;
central_return_air.m_dot = return_air.m_dot - (direct_flow - direct_flow_fresh)*air_density;

direct_return_air = return_air;
direct_return_air.m_dot = (direct_flow - direct_flow_fresh)*air_density;

%flow agregated into air loops for supply side
loop_return_air.m_dot = zone_2_loop*(central_return_air.m_dot.*multiplier);
loop_return_air.w = (zone_2_loop*(central_return_air.m_dot.*central_return_air.w.*multiplier))./loop_return_air.m_dot;
loop_return_air.h = (zone_2_loop*(central_return_air.m_dot.*central_return_air.h.*multiplier))./loop_return_air.m_dot;
%avoid issues with zero or nearly zero flow
A = loop_return_air.m_dot<1e-7 | isnan(loop_return_air.h);
loop_return_air.m_dot(A) = 0;
avg_h = (zone_2_loop*central_return_air.h)./sum(zone_2_loop,2);
loop_return_air.h(A) = avg_h(A);
avg_w = (zone_2_loop*central_return_air.w)./sum(zone_2_loop,2);
loop_return_air.w(A) = avg_w(A);
loop_return_air.T =  psychometric(loop_return_air,'Tdb');
%%heat exchanger with relief flow and fresh air?
end%Ends function return_air_calc

function mixed_air = mix_air(return_air,fresh_air,fresh_flow,supply_flow)
air_density = 1.204; %kg/m^3 %standard air density (1.204 kg/m3) adjusted for the local barometric pressure (standard barometric pressure corrected for altitude, ASHRAE 1997 HOF pg. 6.1).
recycle_air = return_air;
recycle_air.m_dot = (supply_flow - fresh_flow)*air_density;
mixed_air.m_dot = recycle_air.m_dot + fresh_flow*air_density;%mass flow (kg/s) :
mixed_air.w = (recycle_air.m_dot.*recycle_air.w + fresh_flow*air_density.*fresh_air.w)./mixed_air.m_dot;
mixed_air.h = (recycle_air.m_dot.*recycle_air.h + fresh_flow*air_density.*fresh_air.h)./mixed_air.m_dot;
mixed_air.T = psychometric(mixed_air,'Tdb');

% %split off relief flow
% relief_air = return_air;
% relief_air.m_dot = return_air.m_dot - recycle_air.m_dot;
end%Ends function mix_air

function [mixed_air,fresh_air_loop] = economized_outdoor_air(building,min_flow,r_air,fresh_air,loop_flow,l,fan_avail,fan_name,des_day)
fresh_air_loop = min_flow;
T_supply_max = building.hvac.design.Tsupply_h(l);
Ta = fresh_air.T;
Tb = r_air.T;
%% run checks to see if economizer is active
econ_on = false;
if econ_on && isempty(des_day) && T_supply_max>Ta && r_air.T>T_supply_max
    error = 1;
    f = f_index(fan_name,building.fans.name);
    max_flow = building.hvac.outdoor_air.max_flow(l);
    while any(abs(error)>1e-3)
        mixed_air = mix_air(r_air,fresh_air,fresh_air_loop,loop_flow);
        cp = (1006 + 1860*mixed_air.w);%Specific heat in kJ/kg
        [~,~,Q_to_air] = fan_calc(building.fans,fan_name,mixed_air,fan_avail(f));
        T_no_treat = mixed_air.T + Q_to_air./(mixed_air.m_dot.*cp);
        error = T_no_treat - T_supply_max;
        if (error>0 && fresh_air_loop>=max_flow) || (error<0 && fresh_air_loop<=min_flow)
            error = 0;
        end
        d_flow = loop_flow.*(1-(Tb+20)./(error+20))./(Ta-Tb).*error;
        d_flow(isnan(d_flow)) = 0;
        fresh_air_loop = max(min_flow,min(max_flow,fresh_air_loop + d_flow));
    end
else
    mixed_air = mix_air(r_air,fresh_air,fresh_air_loop,loop_flow);
end
end%Ends function economized_outdoor_air

function [loop_T_supply,loop_w_supply,mixed_air_loop,air_nodes] = loop_supply_setpoint(building,profile,air_nodes,fresh_air,fresh_air_loop,loop_return_air,loop_flow,central_flow,loop_2_zone,cp_zone,flow_exiting_zone,net_H2O_gain,T_set,w_set,Q_Tset,T_min,T_max,fan_avail,des_day)
air_density = 1.204; %kg/m^3 %standard air density (1.204 kg/m3) adjusted for the local barometric pressure (standard barometric pressure corrected for altitude, ASHRAE 1997 HOF pg. 6.1).
loop_T_supply = loop_return_air.T;
loop_w_supply = loop_return_air.w;
mixed_air_loop = loop_return_air;
for i = 1:1:length(fresh_air_loop)
    term_i = f_index(i,building.hvac.terminals.loop);
    z_i = building.hvac.terminals.zone(term_i);
    on_loop = f_index(i,building.air_supply_equip.loop);
    if sum(central_flow(z_i)) > 0
        air.T = loop_return_air.T(i); air.w = loop_return_air.w(i); air.h = loop_return_air.h(i); air.m_dot = loop_return_air.m_dot(i); 
        k = f_index({'Fan:ConstantVolume';'Fan:VariableVolume';'Fan:OnOff';'Fan:SystemModel';},building.air_supply_equip.type(on_loop));
        if ~isempty(k)
            fan_name = building.air_supply_equip.name{on_loop(k)};
        else
            k = f_index({'AirLoopHVAC:UnitaryHeatCool';},building.air_supply_equip.type(on_loop));
            u = f_index(building.air_supply_equip.name{on_loop(k)},building.unitary_heat_cool.name);
            fan_name = building.unitary_heat_cool.fan_name{u};
        end
        m_num = f_index(building.air_supply_nodes.name{building.air_supply_equip.outlet{on_loop(end)}},building.manager.node);
        m_type = building.manager.type(m_num);
        if any(strcmp('Scheduled',m_type))
            m_shed = building.manager.schedule{m_num(f_index('Scheduled',m_type))};
            for j = 1:1:length(on_loop)
                out = building.air_supply_equip.outlet{on_loop(j)};
                node = building.air_supply_nodes.name{out};
                if strcmp('AirLoopHVAC:OutdoorAirSystem',building.air_supply_equip.type(on_loop(j)))
                    [mixed_air,fresh_air_loop(i)] = economized_outdoor_air(building,fresh_air_loop(i),air,fresh_air,loop_flow(i),i,fan_avail,fan_name,des_day);
                    mixed_air_loop.T(i) = air.T; mixed_air_loop.w(i) = air.w; mixed_air_loop.h(i) = air.h; mixed_air_loop.m_dot(i) = air.m_dot; 
                    air_nodes.supply_T_set(out) = mixed_air.T;
                else
                    air_nodes.supply_T_set(out) = manager_override(building,profile,node,profile.(m_shed),air,fresh_air,fan_avail);
                end
            end
        else
            f = f_index(fan_name,building.fans.name);
            [~,~,Q_to_air] = fan_calc(building.fans,fan_name,air,fan_avail(f));
            Q_f2zone = Q_to_air*loop_2_zone(z_i,i)./building.zones.multiplier(z_i);%fan heat to zone air supply
            air_nodes.supply_T_set(building.air_supply_equip.inlet{on_loop(1)}) = air.T;
            for j = 1:1:length(on_loop)
                %% making sure simulation of these components is in the order they are on the supply branch
                e = f_index(j,building.air_supply_equip.branch_order(on_loop));
                k = on_loop(e);
                in = building.air_supply_equip.inlet{k};
                out = building.air_supply_equip.outlet{k};
                air_nodes.supply_T_set(out) = air_nodes.supply_T_set(in);
                node = building.air_supply_nodes.name{out};
                switch building.air_supply_equip.type{k}
                    case 'AirLoopHVAC:OutdoorAirSystem'
                        [air,fresh_air_loop(i)] = economized_outdoor_air(building,fresh_air_loop(i),air,fresh_air,loop_flow(i),i,fan_avail,fan_name,des_day);
                        air_nodes.supply_T_set(out) = air.T;
                        mixed_air_loop.T(i) = air.T; mixed_air_loop.w(i) = air.w; mixed_air_loop.h(i) = air.h; mixed_air_loop.m_dot(i) = air.m_dot; 
                    case {'CoilSystem:Cooling:DX';'Coil:Cooling:Water:DetailedGeometry';'Coil:Cooling:Water';'DX:TwoSpeed';'DX:SingleSpeed';}
                        if any(central_flow(z_i).*cp_zone(z_i).*(air.T-T_set.cool(z_i)')>(Q_Tset(z_i,1)-Q_f2zone))%HVAC in cooling mode
                            z_sup_T = Q_Tset(z_i,1)./(cp_zone(z_i).*central_flow(z_i)) + T_set.cool(z_i)';
                            air.T = max(min(z_sup_T(central_flow(z_i)>0)),building.hvac.design.Tsupply_c(i));%constrain supply temperature
                            air_nodes.supply_T_set(out) = manager_override(building,profile,node,air.T,air,fresh_air,fan_avail);
                        end
                        if min(w_set.dehumidify(z_i))<1
                            w_target = psychometric('Tdb',T_set.cool(z_i),'RH',w_set.dehumidify(z_i)*100,'w');
                            w_supply = max(.001,sum((flow_exiting_zone(z_i).*w_target - net_H2O_gain(z_i)))/sum(central_flow(z_i)));% find supply humidity
                            air.w = min(w_supply);
                            if mixed_air_loop.w(i) > air.w
                                air.T =  min(psychometric('Tdb',T_set.cool(z_i),'w',w_supply,'Tdp'));
                                air_nodes.supply_T_set(out) = manager_override(building,profile,node,air.T,air,fresh_air,fan_avail);
                            end
                        end
                    case {'Coil:Heating:Fuel';'Coil:Heating:Electric';'Coil:Heating:Water';}
                        if any(central_flow(z_i).*cp_zone(z_i).*(air.T-T_set.heat(z_i)')<(Q_Tset(z_i,2)-Q_f2zone))%HVAC in heating mode
                            z_sup_T = Q_Tset(z_i,2)./(cp_zone(z_i).*central_flow(z_i)) + T_set.heat(z_i)';
                            T_sup_Theat = min(max(z_sup_T,T_min(z_i)),T_max(z_i));
                            air.T = max(T_sup_Theat(central_flow(z_i)>0));
                            air.T = min(air.T,building.hvac.design.Tsupply_h(i));%constrain supply temperature
                            air_nodes.supply_T_set(out) = manager_override(building,profile,node,air.T,air,fresh_air,fan_avail);
                        end
                    case {'Fan:ConstantVolume';'Fan:VariableVolume';'Fan:OnOff';'Fan:SystemModel'}
                        [~,~,Q_to_air] = fan_calc(building.fans,fan_name,air,fan_avail(f));
                        air.T = air_nodes.supply_T_set(in) + Q_to_air./(air.m_dot.*mean(cp_zone(z_i))/air_density);
                        air_nodes.supply_T_set(out) = air.T;%manager_override(building,profile,node,air.T,air,fresh_air,fan_avail);
                        Q_f2zone = 0;%zero out for any components downstream of fan
                    case 'AirLoopHVAC:UnitaryHeatCool'
                        if any(central_flow(z_i).*cp_zone(z_i).*(air.T-T_set.cool(z_i)')>(Q_Tset(z_i,1)-Q_f2zone))%HVAC in cooling mode
                            z_sup_T = Q_Tset(z_i,1)./(cp_zone(z_i).*central_flow(z_i)) + T_set.cool(z_i)';
                            air.T = min(z_sup_T(central_flow(z_i)>0));
                        elseif any(central_flow(z_i).*cp_zone(z_i).*(air.T-T_set.heat(z_i)')<(Q_Tset(z_i,2)-Q_f2zone))%HVAC in heating mode
                            z_sup_T = Q_Tset(z_i,2)./(cp_zone(z_i).*central_flow(z_i)) + T_set.heat(z_i)';
                            T_sup_Theat = min(max(z_sup_T(central_flow(z_i)>0),T_min(z_i)),T_max(z_i));
                            air.T = min(building.hvac.design.Tsupply_h(i),max(T_sup_Theat));
                        else
                            air = fan_calc(building.fans,fan_name,air,fan_avail(f));
                        end
                        air_nodes.supply_T_set(out) = manager_override(building,profile,building.air_supply_nodes.name{in},air.T,air,fresh_air,fan_avail);
                    case 'Humidifier:Steam:Electric'
                        if max(w_set.humidify(z_i))>0
                            %%need to figure out temperature when humidifying?
                            w_target = psychometric('Tdb',T_set.heat(z_i),'RH',w_set.humidify(z_i)*100,'w');%heating = maybe humidifying
                            w_supply = max(.001,sum((flow_exiting_zone(z_i).*w_target - net_H2O_gain(z_i)))/sum(central_flow(z_i)));% find supply humidity
                            if mixed_air_loop.w(i)<w_supply 
                                air.w = max(w_supply);
                            end
                        end
                        air_nodes.supply_T_set(out) = air.T;
    %                     air_nodes.supply_w_set(out) = manager_override(building,profile,node,air.w,air,fresh_air,fan_avail);
                    otherwise
                        disp('something else on air loop?')
                end
                
            end
        end
        loop_T_supply(i) = air_nodes.supply_T_set(out);
        loop_w_supply(i) = air.w;
    else
        for j = 1:1:length(on_loop)
            out = building.air_supply_equip.outlet{on_loop(j)};
            air_nodes.supply_T_set(out) = loop_return_air.T(i);
        end
    end
end
end%Ends function loop_supply_setpoint

function [T_terminal_supply,unitary_supply_T,mixed_air_zone] = ...
    zones_supply_setpoint(building,central_flow,cp_zone,loop_T_supply,direct_return_air,fresh_air,direct_fresh_air,direct_flow,T_zone_est,T_set,Q_Tset,T_min,T_max,fan_avail)
air_density = 1.204; %kg/m^3 %standard air density (1.204 kg/m3) adjusted for the local barometric pressure (standard barometric pressure corrected for altitude, ASHRAE 1997 HOF pg. 6.1).
n_z = length(building.zones.name);
T_terminal_supply = zeros(n_z,1);
Q_loop_air = zeros(n_z,1);
for z = 1:1:length(building.zones.name)
    term_z = f_index(z,building.hvac.terminals.zone);
    if ~isempty(term_z)
        T_terminal_supply(z) = loop_T_supply(f_index(1,building.hvac.zone_2_loop(:,z)));
        switch building.hvac.terminals.type{term_z}%what if terminals on same loop are of different variety? Not an issue with default buildings
            case 'Uncontrolled'
                %%do nothing
            case 'ConstantVolume:Reheat'
                %%need example
            case 'VAV:Reheat'
                if central_flow(z)>0 && central_flow(z)*cp_zone(z)*(T_terminal_supply(z) - T_set.heat(z))<Q_Tset(z,2)%needs heat
                    c  = f_index(building.hvac.terminals.reheat_coil_name{term_z},building.coils.heating.name);
                    m = building.zones.multiplier(z);
                    Q_reheat_max = min(m*central_flow(z)*cp_zone(z)*(T_max(z) - T_terminal_supply(z)),building.coils.heating.capacity(c));
                    Q_reheat_for_target = m*(Q_Tset(z,2) - central_flow(z)*cp_zone(z)*(T_terminal_supply(z) - T_set.heat(z)));
                    Q_terminal = max(0,min(Q_reheat_for_target,Q_reheat_max));
                    T_terminal_supply(z) = T_terminal_supply(z) + Q_terminal/(central_flow(z)*m*cp_zone(z));
                end
            otherwise
                disp(strcat('air terminal of type: ',building.hvac.terminals.type{term_z}))
        end
        Q_loop_air(z) = central_flow(z)*cp_zone(z)*(T_terminal_supply(z) - T_zone_est(z));
    else
        %not on an air loop (no terminal) check unitary systems
        T_terminal_supply(z) = direct_return_air.T(z);
        Q_loop_air(z) = 0;
    end
end
%% figure out direct flow if unitary system in parallel to AirLoop
Q_direct = Q_Tset+Q_loop_air;
unitary_supply_T = direct_return_air.T;
mixed_air_zone = direct_return_air;
mixed_air_zone.flow = direct_return_air.m_dot/air_density;
for i = 1:1:length(building.unitary_sys.name)
    z = building.unitary_sys.zone(i);
    if direct_flow(z)>0
        r_air.T = direct_return_air.T(z); r_air.w = direct_return_air.w(z); r_air.h = direct_return_air.h(z); r_air.m_dot = direct_return_air.m_dot(z); 
        mixed_air = mix_air(r_air,fresh_air,direct_fresh_air(z),direct_flow(z));%fresh air directly to zone, not outdoor air mixing on loop
        mixed_air_zone.T(z) = mixed_air.T; mixed_air_zone.w(z) = mixed_air.w; mixed_air_zone.h(z) = mixed_air.h; mixed_air_zone.m_dot(z) = mixed_air.m_dot; mixed_air_zone.flow(z) = mixed_air.m_dot/air_density;
        cp_zone = (1006 + 1860*mixed_air.w)*air_density;%Specific heat in kJ/m^3
        fan_name = building.unitary_sys.fan_name{i};
        f = f_index(fan_name,building.fans.name);
        [~,~,Q_to_air] = fan_calc(building.fans,fan_name,mixed_air,fan_avail(f));
        no_treat_zone_air_T = mixed_air.T + Q_to_air./(mixed_air.m_dot.*cp_zone);
        switch building.unitary_sys.type{i}
            case {'PackagedTerminalAirConditioner';'FourPipeFanCoil';}
                if direct_flow(z)*cp_zone*(no_treat_zone_air_T - T_set.cool(z))>Q_direct(z,1)%unitary system in cooling mode
                    T_unit_cooler = min(max(Q_direct(z,1)/(cp_zone*direct_flow(z)) + T_set.cool(z),T_min(z)),T_max(z));
                    unitary_supply_T(z) = min(no_treat_zone_air_T,T_unit_cooler);
                elseif direct_flow(z)*cp_zone*(no_treat_zone_air_T - T_set.heat(z))<Q_direct(z,2)%unitary system  in heating mode
                    T_unit_heater = min(Q_direct(z,2)/(cp_zone*direct_flow(z)) + T_set.heat(z),T_max(z));
                    unitary_supply_T(z) = max(no_treat_zone_air_T,T_unit_heater);
                end
            case 'UnitHeater'
                if direct_flow(z)*cp_zone*(no_treat_zone_air_T - T_set.heat(z))<Q_direct(z,2)%unitary system  in heating mode
                    T_unit_heater = min(Q_direct(z,2)/(cp_zone*direct_flow(z)) + T_set.heat(z),T_max(z));
                    unitary_supply_T(z) = max(no_treat_zone_air_T,T_unit_heater);
                end
            case 'WindowAirConditioner'
                disp('need WindowAirConditioner')
            otherwise
                disp('need other types of unitary system')
        end
    end
end
end%Ends function zones_on_loop

function fresh_air_zone = design_zone_outdoor_air(building,occupancy)
%% Determine max fresh air requirements
fresh_air_zone = zeros(length(building.zones.name),1);
for i = 1:1:length(building.setpoints.name)
    z = building.setpoints.zone(i);
    switch building.setpoints.outdoor_flow_method{i}
        case 'Flow/Person'
            fresh_air_zone(z) = building.setpoints.outdoor_flow_value(i)*occupancy(z);
        case 'Flow/Area'
            fresh_air_zone(z) = building.setpoints.outdoor_flow_value(i)*building.zones.floor_area(z);
        case 'Flow/Zone'
            fresh_air_zone(z) = building.setpoints.outdoor_flow_value(i);
    end
end
end%Ends function design_zone_outdoor_air

function [central_flow,direct_flow,supply_flow,direct_fresh_air,fresh_air_loop] = des_day_zone_flow(building,fresh_air_zone,cp_supply,Q_Tset,T_set,T_min,T_max,flow_imbalance)
zone_flow_req = zeros(length(T_min),1);
h = Q_Tset(:,2)>0 & ~T_set.no_hvac';
zone_flow_req(h) = Q_Tset(h,2)./(cp_supply(h).*(T_max(h) - T_set.heat(h)'));
c =  Q_Tset(:,1)<0 & ~T_set.no_hvac';
zone_flow_req(c) =  Q_Tset(c,1)./(cp_supply(c).*(T_min(c) - T_set.cool(c)'));
zone_flow_req = max(zone_flow_req, -flow_imbalance.*(~T_set.no_hvac'));%volumetric flow m^3/s

zone_flow = max(fresh_air_zone,zone_flow_req);%volumetric flow m^3/s
direct_flow = zeros(length(building.zones.name),1);%flow for zones not on a HVAC loop (unit heaters)
for i = 1:1:length(building.unitary_sys.name)
    l = f_index(building.unitary_sys.outlet{i},building.hvac.loop.supply_outlet);
    if isempty(l)
        z = building.unitary_sys.zone(i);
        direct_flow(z) = zone_flow(z);
    end
end
central_flow = zone_flow - direct_flow;
supply_flow = building.hvac.zone_2_loop*(central_flow.*building.zones.multiplier);%mass flow for the loop (kg/s)
fresh_air_loop = building.hvac.zone_2_loop*fresh_air_zone;
direct_fresh_air = fresh_air_zone - fresh_air_loop*loop_2_zone;
end%Ends function des_day_zone_flow

function [T_min,T_max] = zone_limits(building,T_zone)
T_min = zeros(length(T_zone),1);
T_max = zeros(length(T_zone),1);
for i = 1:1:length(building.setpoints.name)
    z = building.setpoints.zone;
    switch building.setpoints.cooling_T_method{i}
        case 'SupplyAirTemperature'
            T_min(z) = building.setpoints.cooling_T_des(i);
        otherwise
            T_min(z) = T_zone(z) - building.setpoints.cooling_dT_des(i);
    end
    switch building.setpoints.heating_T_method{i}
        case 'SupplyAirTemperature'
            T_max(z) = building.setpoints.heating_T_des(i);
        otherwise
            T_max(z) = T_zone(z) + building.setpoints.heating_dT_des(i);
    end
end
for i = 1:1:length(building.manager.name)
    switch building.manager.type{i}
        case 'SingleZone:Reheat'
            z = f_index(building.manager.zone{i},building.zones.name);
            T_min(z) = building.manager.min_temp(i);
            T_max(z) = building.manager.max_temp(i);
    end
end
end%Ends function zone_limits