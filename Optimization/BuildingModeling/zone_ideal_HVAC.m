function [Q_zone,zone_flow,T_supply,T_loop_supply,m_v_supply,fresh_air_zone,zone_exhaust_fan,loop_mixed_air,non_loop_mixed_air,Q_terminal] = ...
    zone_ideal_HVAC(building,profile,mixing,infiltration,Q_zone_Tcool,Q_zone_Theat,T_set,T_zone_est,m_v_zone,T_air,m_v_air,occupancy,des_day)
%%Determine water vapor mass fraction and specific heats (need to update)
air_density = 1.225; %kg/m^3
n = length(building.zones.name);
m_v_supply = m_v_air*ones(n,1);
cp_supply = 1006 + 1860*m_v_supply;

loops = length(building.HVAC.zone_2_loop(:,1));
multiplier = building.zones.multiplier;

%%exhaust fan
n_f = length(building.HVAC.fans.name);
zone_exhaust_fan = zeros(1,n);
exhaust_flow = zeros(n,1);
for i = 1:1:n
    k = nonzeros((1:n_f)'.*strcmpi(building.zones.exhaust{i},building.HVAC.fans.inlet));
    if ~isempty(k)
        sched = building.HVAC.fans.schedule{k};
        if ~isnan(building.HVAC.fans.flow_rate(k))
            exhaust_flow(i) = building.HVAC.fans.flow_rate(k)*profile.(sched);%volumetric flow (m^3/s)  %See page 1601 of the reference manual
            zone_exhaust_fan(i)= exhaust_flow(i).*building.HVAC.fans.pressure_rise(k)/(building.HVAC.fans.fan_efficiency(k)*air_density);% kg/s * Pa / (kg/m^3) = N*m/s = W
        end
    end
    exhaust_flow(i) = exhaust_flow(i) - sum(mixing(:,i));%balanced exhaust flow
end
%%Determine minimum mass flow when actively heating/cooling the zone
%When heating, provide exactly Q_sys and at least fresh air flow. Either m_sys>fresh air & T_supply = T_set, or m_sys = fresh air & T_supply<T_set_heat
%%logic of heating/cooling/passive ventilation
zone_mode = cell(length(T_zone_est),1);
loop_mode = cell(loops,1);
zone_mode(Q_zone_Theat>0 & ~T_set.no_hvac') = {'heat'};
zone_mode(Q_zone_Tcool<0 & ~T_set.no_hvac') = {'cool'};
zone_mode((Q_zone_Tcool>0 & Q_zone_Theat<0) | T_set.no_hvac') = {'passive'};
Q_zone = Q_zone_Theat;%default is the lower zone temperature this way a passive system will not overcool a zone (reheat coils)
min_zone_flow = zeros(n,1);
for i = 1:1:n
    switch zone_mode{i}
        case 'heat'
            min_zone_flow(i) = Q_zone_Theat(i)/(cp_supply(i)*(building.zones.t_supply_h(i) - T_set.heat(i)));
        case 'cool'
            min_zone_flow(i) = Q_zone_Tcool(i)/(cp_supply(i)*(building.zones.t_supply_c(i) - T_set.cool(i)));
            Q_zone(i) = Q_zone_Tcool(i);
    end  
end
min_zone_flow = max(min_zone_flow,-air_density*(infiltration - exhaust_flow - sum(mixing,1)' + sum(mixing,2)).*(~T_set.no_hvac'));%mass flow kg/s


T_zone_est2 = T_zone_est;
T_loop_supply = zeros(loops,1);
m = (Q_zone_Tcool - Q_zone_Theat)./(T_set.cool - T_set.heat)';
for repeat = 1:1:2
    return_air_difference = building.HVAC.zone_2_loop*(air_density*(infiltration - exhaust_flow - sum(mixing,1)' + sum(mixing,2)).*multiplier);%mass flow that must be made up by fresh air(kg/s)
    [zone_flow,fresh_air_zone,supply_flow,fresh_air_loop] = air_mass_balance(building,profile,min_zone_flow,zone_mode,occupancy,return_air_difference,Q_zone,T_zone_est2,m_v_zone,multiplier,des_day);
    return_air_flow = zone_flow + air_density*(infiltration - exhaust_flow - sum(mixing,1)' + sum(mixing,2));%mass flow (kg/s)
    loop_mixed_air = outdoor_air_mixing(supply_flow,fresh_air_loop,return_air_flow,building.HVAC.zone_2_loop,building.zones.multiplier,T_zone_est,m_v_zone,T_air,m_v_air);%% determine mixed air flow temperature and humidity
    m_v_supply = building.HVAC.zone_2_loop'*loop_mixed_air.w;
    cp_supply = 1006 + 1860*m_v_supply;

    for i = 1:1:loops
        term_i = nonzeros((1:length(building.HVAC.terminals.name))'.*(building.HVAC.terminals.loop == i));
        z_i = building.HVAC.terminals.zone(term_i);
        T_min_treatment = max(min(loop_mixed_air.T(i),building.HVAC.design.Tsupply_h(i)),building.HVAC.design.Tsupply_c(i));
        if any(zone_flow(z_i).*cp_supply(z_i).*(T_min_treatment-T_set.cool(z_i)')>Q_zone_Tcool(z_i))%HVAC in cooling mode
            loop_mode(i) = {'cool'};
            T_sup_Tcool = min(max(Q_zone_Tcool(z_i)./(cp_supply(z_i).*zone_flow(z_i)) + T_set.cool(z_i)',building.zones.t_supply_c(z_i)),building.zones.t_supply_h(z_i));
            T_loop_supply(i) = min(T_sup_Tcool);%%for multiple zones on same loop (VAV boxes)
% %                 T_loop_supply(i) = max(building.HVAC.design.Tsupply_c(i),min(T_sup_Tcool));
        elseif any(zone_flow(z_i).*cp_supply(z_i).*(T_min_treatment-T_set.heat(z_i)')<Q_zone_Theat(z_i))%HVAC in heating mode
            loop_mode(i) = {'heat'};
            T_sup_Theat = min(max(Q_zone_Theat(z_i)./(cp_supply(z_i).*zone_flow(z_i)) + T_set.heat(z_i)',building.zones.t_supply_c(z_i)),building.zones.t_supply_h(z_i));
            T_loop_supply(i) = min(building.HVAC.design.Tsupply_h(i),max(T_sup_Theat));
        else
            loop_mode(i) = {'passive'};
            T_loop_supply(i) = T_min_treatment;
            for j = 1:1:length(building.HVAC.loop.component_name{i})
                if any(strcmpi(building.HVAC.loop.component_type{i}{j},{'Fan:ConstantVolume';'Fan:VariableVolume';'Fan:OnOff'}))
                    air_est.T = T_min_treatment; air_est.w = loop_mixed_air.w(i);air_est.m_dot = loop_mixed_air.m_dot(i);air_est.h = psychometric_h(air_est);                    
                    f = nonzeros((1:length(building.HVAC.fans.name))'.*strcmp(building.HVAC.loop.component_name{i}{j},building.HVAC.fans.name));
                    air_out = fan_calc(building.HVAC.fans,building.HVAC.loop.component_name{i}{j},'operation',air_est,profile.(building.HVAC.fans.schedule{f}));
                    T_loop_supply(i) = air_out.T;
                end
            end
        end
    end
    T_split_supply = building.HVAC.zone_2_loop'*T_loop_supply;
    [T_split_supply,m_v_supply,non_loop_mixed_air] = unitary_systems(building.HVAC.unitary_sys,Q_zone_Theat,T_split_supply,m_v_supply,zone_flow,T_zone_est,m_v_zone,fresh_air_zone,zone_mode,T_air,m_v_air,T_set,building.zones.t_supply_h);%Unit heaters/coolers in zones not part of central air loops
    
    Q_terminal = max(0,Q_zone_Theat-zone_flow.*cp_supply.*(T_split_supply-T_set.heat')) - max(0,zone_flow.*cp_supply.*(T_split_supply-T_set.cool') - Q_zone_Tcool);%If reheat is needed to stay above T_set.heat
    T_zone_est2 = (m.*T_set.heat' + zone_flow.*cp_supply.*T_split_supply - Q_zone_Theat + Q_terminal)./(m + zone_flow.*cp_supply);%estimate of zone temperature with air treatment
    T_zone_est2(T_set.cool == T_set.heat) = T_set.cool(T_set.cool == T_set.heat);
    Q_zone = zone_flow.*cp_supply.*(T_split_supply - T_zone_est2) + Q_terminal;%zone energy exchange   
end
T_supply = Q_zone./(cp_supply.*zone_flow) + T_zone_est2;
T_supply(zone_flow == 0) = 0;
Q_zone(zone_flow == 0) = 0;
end%Ends function zone_ideal_HVAC


function [zone_flow,fresh_air,supply_flow,fresh_air_loop] = air_mass_balance(building,profile,zone_flow,zone_mode,occupancy,return_air_difference,Q_zone,T_zone_est,m_v_zone,multiplier,des_day)
%%air mass balance
%Determine flow rate back to zones, based on fan flow limits, fresh air requirements and if HVAC is on
loops = length(building.HVAC.zone_2_loop(:,1));
air_density = 1.225; %kg/m^3
cp_zone = 1006 + 1860*m_v_zone; %J/kg*K
m_zone = zeros(loops,1);
min_cooling_T = zeros(loops,1);
unit_flow = 0*zone_flow;%flow for zones not on a HVAC loop (unit heaters)
fresh_air_zone = 0*zone_flow;
if isempty(des_day)%%do VAV terminal air calculations   
    term_flow = zeros(length(building.HVAC.terminals.name),1); 
    for i = 1:1:loops
        term_i = nonzeros((1:length(building.HVAC.terminals.name))'.*(building.HVAC.terminals.loop == i));
        %Cooling mode:supply temperature is based on limiting zone, other zones operate with reduced flow. If a terminal minimum flow is hit, re-heat prevents from overcooling
        %Heating mode: minimum flow. re_heat does all necessary heating
        z_i = building.HVAC.terminals.zone(term_i);
        [min_cooling_T(i),m_zone(i)] = min(max(building.zones.t_supply_c(z_i),T_zone_est(z_i) + Q_zone(z_i)./(cp_zone(z_i).*building.HVAC.terminals.max_flow(term_i))));
        for j = 1:1:length(term_i)
            if strcmp(building.HVAC.terminals.type,'Uncontrolled')
                term_flow(term_i(j)) = zone_flow(z_i(j));
            elseif strcmp(zone_mode{term_i(j)},'cool')
                if j == m_zone(i)
                    term_flow(term_i(j)) = building.HVAC.terminals.max_flow(term_i(j));
                else
                    term_flow(term_i(j)) = max(building.HVAC.terminals.min_flow(term_i(j)),Q_zone(z_i(j))./(cp_zone(z_i(j)).*(min_cooling_T(i) - T_zone_est(z_i(j)))));
                end
            elseif strcmp(zone_mode{term_i(j)},'heat')
                term_flow(term_i(j)) = building.HVAC.terminals.min_flow(term_i(j));
            elseif strcmp(zone_mode{term_i(j)},'passive')
                term_flow(term_i(j)) = building.HVAC.terminals.min_flow(term_i(j));
            end
        end
    end
    zone_flow(building.HVAC.terminals.zone) = term_flow;
    %%unitary systems not on HVAC loop
    if ~isempty(building.HVAC.unitary_sys)
        for i = 1:1:length(building.HVAC.unitary_sys.name)
            if building.HVAC.unitary_sys.loop(i) == 0
                z = building.HVAC.unitary_sys.zone(i);
                f = nonzeros((1:length(building.HVAC.fans.name))'.*strcmp(building.HVAC.unitary_sys.fan_name{i},building.HVAC.fans.name));
                unit_flow(z) = building.HVAC.fans.flow_rate(f)*profile.(building.HVAC.fans.schedule{i});
            end
        end
        
        for i = 1:1:length(building.HVAC.outdoor_air.name)%loop through outside_air controllers to get minimum
            z = building.HVAC.outdoor_air.zone(i);
            if z>0
                fresh_air_zone(z) = outdoor_air.min_flow(i)*profile.(building.HVAC.outdoor_air.schedule{i});
            end
        end
    end
    
else
    %% Determine max fresh air requirements
    % If outside air method is flow/zone, the input outside air flow per zone value will be used, even
    % if it is zero or blank. If outside air method is sum, the sum of the outside air flow per person *
    % DesignNumberOfPeople + outside air flow per area * ZoneArea will be used. If outside air method
    % is maximum, the maximum of the outside air flow per person * DesignNumberOfPeople and outside
    % air flow per area * ZoneArea will be used. If outside air method is flow/person, outside air flow per
    % person will be used to calculate the design minimum outside airflow rate.
    fresh_air_zone = building.zones.outdoor_flow.value*air_density;%mass flow (kg/s
    A = strcmp('Flow/Person',building.zones.outdoor_flow.method);
    fresh_air_zone(A) = building.zones.outdoor_flow.value(A).*occupancy(A)*air_density;%mass flow (kg/s
    B = strcmp('Flow/Area',building.zones.outdoor_flow.method);
    fresh_air_zone(B) = building.zones.outdoor_flow.value(B).*building.zones.floor_area(B)*air_density;%mass flow (kg/s
    zone_flow = max(fresh_air_zone,zone_flow);
    if ~isempty(building.HVAC.unitary_sys)
        for i = 1:1:length(building.HVAC.unitary_sys.name)
            if building.HVAC.unitary_sys.loop(i) == 0
                unit_flow(building.HVAC.unitary_sys.zone(i)) = zone_flow(building.HVAC.unitary_sys.zone(i));
            end
        end
    end
end
supply_flow = building.HVAC.zone_2_loop*(zone_flow.*multiplier);%mass flow (kg/s)
fresh_air_loop = max(0,-return_air_difference);
if isempty(des_day)
    [minimum_flow,limit_flow,fresh_air_loop] = supply_flow_limit(building.HVAC,profile,building.HVAC.outdoor_air,fresh_air_loop,loops);%enforces minimum and maximum flow rates
    supply_flow = max(minimum_flow*air_density,min(limit_flow*air_density,supply_flow));
end
loop_2_zone = 0*building.HVAC.zone_2_loop';
for i = 1:1:loops
    if supply_flow(i)>0 %split flow based on heating/cooling request flow (avoid NaN from dividing by zero)
        min_loop_flow = (zone_flow'*(building.HVAC.zone_2_loop(i,:)'.*multiplier));
        if min_loop_flow>0
            loop_2_zone(:,i) = (zone_flow.*building.HVAC.zone_2_loop(i,:)'.*multiplier)/min_loop_flow;
        else %avoid NaN, split to zones based on max zone flow
            zone_request = 0*zone_flow;
            z = nonzeros((1:length(zone_flow)).*building.HVAC.zone_2_loop(i,:));
            zone_request(z) = building.HVAC.terminals.max_flow(building.HVAC.terminals.zone == z).*multiplier(z)*air_density;%mass flow (kg/s)
            loop_2_zone(:,i) = zone_request/sum(zone_request);
        end
    end
end
zone_flow = loop_2_zone*supply_flow./multiplier + unit_flow;%mass flow (kg/s)
fresh_air = loop_2_zone*fresh_air_loop./multiplier + fresh_air_zone;%mass flow (kg/s)
end%Ends function air_mass_balance

function [minimum_flow,limit_flow,fresh_air_loop] = supply_flow_limit(hvac,profile,outdoor_air,fresh_air_loop,loops)
air_density = 1.225; %kg/m^3
for i = 1:1:length(outdoor_air.name)
    l = outdoor_air.loop(i);
    if l>0
        fresh_air_loop(i) = max(fresh_air_loop(i),air_density*outdoor_air.min_flow(i)*profile.(outdoor_air.schedule{i}));%mass flow (kg/s)
    end
end
    
%reduce flow if air_terminal cannot meet supply request

%reduce flow if any equipment in air_loop limits flow
limit_flow = air_density*hvac.design.max_loop_flow;%mass flow (kg/s)
minimum_flow = fresh_air_loop;%mass flow (kg/s)
%%go through equipment on loop and reduce if schedule is reduced
for i = 1:1:loops
    for j = 1:1:length(hvac.loop.component_name{i})
        type = hvac.loop.component_type{i}{j};
        name = hvac.loop.component_name{i}{j};
        switch type
             case {'CoilSystem:Cooling:DX';'Coil:Cooling:Water:DetailedGeometry';'Coil:Cooling:Water'}
                if ~any(strcmpi(name,hvac.coils.cooling.name))
                    name = hvac.coil_system.coil_name{nonzeros((1:length(hvac.coil_system.name))'.*strcmpi(name,hvac.coil_system.name))};
                end
                k2 = nonzeros((1:length(hvac.coils.cooling.name))'.*strcmpi(name,hvac.coils.cooling.name));
                sched = hvac.coils.cooling.schedule{k2};
                limit_flow(i) = min(limit_flow(i),air_density*hvac.coils.cooling.rated_air_flow(k2)*profile.(sched));%mass flow (kg/s)
            case {'Fan:ConstantVolume';'Fan:OnOff';}
                k2 = nonzeros((1:length(hvac.fans.name))'.*strcmpi(name,hvac.fans.name));
                sched = hvac.fans.schedule{k2};
                limit_flow(i) = air_density*hvac.fans.flow_rate(k2)*profile.(sched);%mass flow (kg/s)
                minimum_flow(i) =air_density*hvac.fans.flow_rate(k2)*profile.(sched);%mass flow (kg/s)
            case 'Fan:VariableVolume'
                k2 = nonzeros((1:length(hvac.fans.name))'.*strcmpi(name,hvac.fans.name));
                sched = hvac.fans.schedule{k2};
                limit_flow(i) = min(limit_flow(i),air_density*hvac.fans.flow_rate(k2)*profile.(sched));%mass flow (kg/s)
                minimum_flow(i) = max(minimum_flow(i),air_density*hvac.fans.min_flow_rate(k2)*profile.(sched));%mass flow (kg/s)                
        end     
    end
end
end%Ends function supply_flow_limit

function mixed_air = outdoor_air_mixing(supply_flow,fresh_air_loop,return_air_flow,zone_2_loop,multiplier,T_zone,m_v_zone,T_air,m_v_air)
%note: should upgrade to include economizer (none in initial default buildings)
%%The outside air mixer splits the primary air system (AirLoopHVAC) return air into relief and recirculated air streams.
%%Then it mixes the outside air stream with the recirculated air stream to obtain the mixed air stream

%determine air returning to HVAC from zone, and agregate into loops
return_air_zone.m_dot = return_air_flow;
return_air_zone.T = T_zone;
return_air_zone.w = ones(length(T_zone),1).*m_v_zone;
return_air_zone.h = psychometric_h(return_air_zone);

%flow agregated into air loops for supply side
return_loop_air.m_dot = zone_2_loop*(return_air_zone.m_dot.*multiplier);
return_loop_air.w = (zone_2_loop*(return_air_zone.m_dot.*return_air_zone.w.*multiplier))./return_loop_air.m_dot;
return_loop_air.h = (zone_2_loop*(return_air_zone.m_dot.*return_air_zone.h.*multiplier))./return_loop_air.m_dot;
%avoid issues with zero or nearly zero flow
A = return_loop_air.m_dot<1e-7 | isnan(return_loop_air.h);
return_loop_air.m_dot(A) = 0;
avg_h = (zone_2_loop*return_air_zone.h)./sum(zone_2_loop,2);
return_loop_air.h(A) = avg_h(A);
avg_w = (zone_2_loop*return_air_zone.w)./sum(zone_2_loop,2);
return_loop_air.w(A) = avg_w(A);
return_loop_air.T =  psychometric_T(return_loop_air);

%split off relief flow
relief_air = return_loop_air;
relief_air.m_dot = return_loop_air.m_dot - (supply_flow-fresh_air_loop);
return_loop_air.m_dot = supply_flow-fresh_air_loop;

fresh_air.m_dot = fresh_air_loop;%mass flow (kg/s)
fresh_air.T = T_air;
fresh_air.w = m_v_air;
fresh_air.h = psychometric_h(fresh_air);

%%heat exchanger with relief flow and fresh air?
    
%%

mixed_air.m_dot = supply_flow;%mass flow (kg/s) :
mixed_air.h = (return_loop_air.m_dot.*return_loop_air.h + fresh_air.m_dot.*fresh_air.h)./mixed_air.m_dot;
mixed_air.w = (return_loop_air.m_dot.*return_loop_air.w + fresh_air.m_dot.*fresh_air.w)./mixed_air.m_dot;
%avoid divide by zeros
A = mixed_air.m_dot<1e-8 | isnan(mixed_air.h);
mixed_air.m_dot(mixed_air.m_dot<1e-8) = 0;%mass flow (kg/s)
mixed_air.h(A) = fresh_air.h(1);
mixed_air.w(A) = fresh_air.w(1);
mixed_air.T = psychometric_T(mixed_air);
end%Ends function outdoor_air

function [T_supply,m_v_supply,non_loop_mixed_air,zone_mode] = unitary_systems(unit,Q_zone_Theat,T_supply,m_v_supply,zone_flow,T_zone,m_v_zone,fresh_air_zone,zone_mode,T_air,m_v_air,T_set,T_max)
%% Unitary systems that are in a zone, not on an HVAC loop
non_loop_mixed_air.T = [];
non_loop_mixed_air.w = [];
j = 0;
if ~isempty(unit)
    for i = 1:1:length(unit.name)
        if unit.loop(i) == 0 %not on a loop
            z = unit.zone(i);
            j = j+1;
            mixed_T = T_air;
            m_v_supply(z) = m_v_air;
            T_supply(z) = mixed_T;
            h = psychometric_h(T_air,m_v_air);
            if zone_flow(z)>0
                h = ((zone_flow(z) - fresh_air_zone(z)).*psychometric_h(T_zone(z),m_v_zone(z)) + fresh_air_zone(z).*psychometric_h(T_air,m_v_air))./zone_flow(z);
                m_v_supply(z) = ((zone_flow(z) - fresh_air_zone(z)).*m_v_zone(z) + fresh_air_zone(z).*m_v_air)./zone_flow(z);
                mixed_T = psychometric_T(h,m_v_supply(z));
                cp_supply = 1006 + 1860*m_v_supply(z);
                T_supply(z) = mixed_T;
                if strcmp(unit.type{i},'heat_cool')
                    if zone_flow(z)*cp_supply*(mixed_T - T_set.cool(z))>Q_zone_Tcool(z)%unitary system in cooling mode
                        zone_mode(z) = {'cool'};
                        T_supply(z) = min(max(Q_zone_Tcool(z)/(cp_supply*zone_flow(z)) + T_set.cool(z),building.zones.t_supply_c(z)),building.zones.t_supply_h(z));
                    elseif zone_flow(z)*cp_supply*(mixed_T - T_set.heat(z))<Q_zone_Theat(z)%unitary system  in heating mode
                        zone_mode(z) = {'heat'};
                        T_supply(z) = min(Q_zone_Theat(z)./(cp_supply.*zone_flow(z)) + T_set.heat(z),T_max(z));
                    end
                elseif strcmp(unit.type{i},'heater')
                    if zone_flow(z)*cp_supply*(mixed_T - T_set.heat(z))<Q_zone_Theat(z)%unitary system in heating mode
                        zone_mode(z) = {'heat'};
                        T_supply(z) = min(Q_zone_Theat(z)./(cp_supply.*zone_flow(z)) + T_set.heat(z),T_max(z));
                    end
                end
            end
            non_loop_mixed_air.T(j,1) = mixed_T;
            non_loop_mixed_air.w(j,1) = m_v_supply(z);
            non_loop_mixed_air.m_dot(j,1) = zone_flow(z);
            non_loop_mixed_air.h(j,1) = h;
        end
    end
end
end%Ends function unitary_systems