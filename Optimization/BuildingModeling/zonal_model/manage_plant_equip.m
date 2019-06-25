function [plant_nodes,e_use,T_tank,tank_on] = manage_plant_equip(building,profile,plant_nodes,e_use,dt,T_air,w_air,T_mains,T_tank,tank_on)
%need to set half loop order for any demand dependacies on supply loops 
%%then can eliminate difference between plant and condenser  loops 
for i = 1:1:length(building.plant_loop.name)
    l = f_index(i,building.plant_loop.loop_order);
    T_tank1 = plant_nodes.demand_temperature(building.plant_demand_nodes.inlet_node(l));%half loop tank temperature from previous time step (pg 468)
    T_tank2 = plant_nodes.supply_temperature(building.plant_supply_nodes.inlet_node(l));%half loop tank temperature from previous time step (pg 468)
    %% Plant demand loops
    [plant_nodes,pump_power,pump_heat] = flow_resolver_demand(building,plant_nodes,profile,l,T_mains,T_tank1);
    e_use.pumps = e_use.pumps + pump_power;%Watts
    plant_nodes = demand_loop_energy_balance(building,profile,T_mains,plant_nodes,l);
    n = building.plant_demand_nodes.outlet_node(l);
    loop_flow_rate = plant_nodes.demand_flow(n);
    T_return = half_loop_tank(building.plant_loop.loop_volume(l)/2,building.plant_loop.fluid_density(l),loop_flow_rate,T_tank2,pump_heat,plant_nodes.demand_temperature(n),dt);% connect demand half loop to supply half loop
    plant_nodes.supply_temperature(building.plant_supply_nodes.inlet_node(l)) = T_return;
    %% Plant supply loops
    plant_nodes.supply_flow = flow_resolver_supply(building,plant_nodes.supply_flow,loop_flow_rate,T_return,building.plant_loop.load_scheme{l},profile,l);
    [plant_nodes,e_use,T_tank,tank_on] = supply_loop(building,profile,plant_nodes,e_use,T_tank,tank_on,T_air,w_air,l,dt);
    n = building.plant_supply_nodes.outlet_node(l);
    plant_nodes.demand_temperature(building.plant_demand_nodes.inlet_node(l)) = ...
        half_loop_tank(building.plant_loop.loop_volume(l)/2,building.plant_loop.fluid_density(l),plant_nodes.supply_flow(n),T_tank1,0,plant_nodes.supply_temperature(n),dt);% connect supply half loop to demand half loop  
end
end%Ends function manage_plant_equip

function [node,pump_power,pump_heat] = flow_resolver_demand(building,node,profile,l,T_mains,T_tank1)
%% resolve loop flows
%see details starting on pg 461
d_equip = building.plant_demand_equip;
on_loop = f_index(l,d_equip.loop);
%%connect flow requests from individual components
for i = 1:1:length(on_loop)
    j = on_loop(i);
    in = d_equip.inlet{j};
    out = d_equip.outlet{j};
    switch d_equip.type{j}
        case 'Mixer'
            node.demand_flow(out) = sum(node.demand_flow(in));
        case 'Splitter'
            node.demand_flow(in) = sum(node.demand_flow(out));
        case 'WaterUse:Connections'
            we = f_index(1,strcmp(d_equip.name{j},building.water_use.name) & strcmp('Equipment',building.water_use.type));
            flow = water_equipment(building,profile,T_mains,T_tank1,we);
            node.demand_flow(out) = flow;
            node.demand_flow(in) = node.demand_flow(out);
        otherwise%flow is already set by component (coil)
            node.demand_flow(out) = max(node.demand_flow(in),node.demand_flow(out));%branch flow determined by component with largest request
            node.demand_flow(in) = node.demand_flow(out);
    end
end
loop_flow_request = node.demand_flow(building.plant_demand_nodes.outlet_node(l))/building.plant_loop.fluid_density(l);% convert mass flow to volumetric flow (parameters saved in building file in volumetric flow)

%find pump & supply limits
pump_flow = 0;
pump_heat = 0;
pump_power = 0;
demand_bypass_flow = 0;
if any(strncmp('Pump:',d_equip.type,5) & d_equip.loop == l)%A Pump, if present, in the demand side of plant loop must be the first component of the inlet branch
    k = f_index(1,strncmp('Pump:',d_equip.type,5) & d_equip.loop == l);
    p = f_index(d_equip.name{k},building.pump.name);
else
    k = f_index(1,strncmp('Pump:',building.plant_supply_equip.type,5) & building.plant_supply_equip.loop == l);
    p = f_index(building.plant_supply_equip.name{k},building.pump.name);
    switch building.pump.type{p}
        case 'VariableSpeed'
            if loop_flow_request>0 && loop_flow_request<building.pump.design_min_flow(p)
                demand_bypass_flow = building.pump.design_min_flow(p) - loop_flow_request;
                pump_flow = building.pump.design_min_flow(p);
            elseif loop_flow_request>building.pump.design_flow(p)
                pump_flow = building.pump.design_flow(p);
            else
                pump_flow = loop_flow_request;
            end
            PLR = pump_flow/building.pump.design_flow(p);
            frac_power = building.pump.c1(p) + building.pump.c2(p)*PLR + building.pump.c3(p)*PLR^2 + building.pump.c4(p)*PLR^3;
        case 'ConstantSpeed'
            if ~isempty(building.pump.schedule{p})
                sched = profile.(building.pump.schedule{p});
            elseif loop_flow_request>0
                sched = 1;
            else
                sched = 0;
            end
            frac_power = sched;
            pump_flow = building.pump.design_flow(p)*sched;
            demand_bypass_flow = pump_flow - loop_flow_request;
    end
    pump_power = building.pump.design_power(p)*frac_power;%Watts
    shaft_power = building.pump.motor_eff(p)*pump_power;%Watts
    pump_heat = shaft_power + (pump_power - shaft_power)*building.pump.fluid_ineff_frac(p);
end
%%put correct flows on demand side branches
node.demand_flow(building.plant_demand_nodes.inlet_node(l)) = pump_flow*building.plant_loop.fluid_density(l);
for i = 1:1:length(on_loop)
    j= on_loop(i);
    in = d_equip.inlet{j};
    out = d_equip.outlet{j};
    switch d_equip.type{j}
        case 'Mixer'
            node.demand_flow(out) = sum(node.demand_flow(in));
        case 'Splitter'
            for k = 1:1:length(on_loop)
                if any(d_equip.inlet{on_loop(k)} == out)
                    if d_equip.bypass(on_loop(k))
                        node.demand_flow(d_equip.inlet{on_loop(k)}) = demand_bypass_flow.*building.plant_loop.fluid_density(l);% convert back to mass flow
                    end
                end
            end
        otherwise
            node.demand_flow(out) = node.demand_flow(in);
    end
end

node.supply_flow(building.plant_supply_nodes.inlet_node(l)) = node.demand_flow(building.plant_demand_nodes.outlet_node(l));% connect demand half loop to supply half loop
end%Ends function flow_resolver_demand

function supply_flow = flow_resolver_supply(building,supply_flow,flow,T_return,method,profile,l)
%%supply side
Cp_water = 4186; %J/kg*K
load = flow*Cp_water*(T_return - building.plant_loop.exit_temperature(l));% (Watts) kg/s*J/kgK * K = J/s
on_loop = f_index(l,building.plant_supply_equip.loop);
for i = 1:1:length(on_loop)
    in = building.plant_supply_equip.inlet{on_loop(i)};
    out = building.plant_supply_equip.outlet{on_loop(i)};
    switch building.plant_supply_equip.type{on_loop(i)}
        case 'Mixer'
            supply_flow(out) = sum(supply_flow(in));
        case 'Splitter'
            supply_flow(out) = load_distribution(building,out,supply_flow(in),T_return,load,method,profile);
        case {'Pump:VariableSpeed','Pump:ConstantSpeed'}
            supply_flow(out) = flow;
        otherwise
            supply_flow(out) = supply_flow(in);
    end
end
end%Ends function flow_resolver_supply

function T_new = half_loop_tank(volume,den,flow,T_old,pump,T_inlet,dt)
if flow>0
    Cp_water = 4186; %J/kg*K
    M_tank = volume*den;%Mass of fluid in tank
    A = (flow*Cp_water*T_inlet+pump)/(flow*Cp_water);%Watts / J/s*K
    T_new = (T_old - A)*exp(-flow*Cp_water*dt/(M_tank*Cp_water)) + A;
else
    T_new = T_old;
end
end%Ends function half_loop_tank

function node = demand_loop_energy_balance(building,profile,T_mains,node,l)
Cp_water = 4186; %J/kg*K
equip = building.plant_demand_equip;
on_loop = f_index(l,equip.loop);
for i = 1:1:length(on_loop)
    in = equip.inlet{on_loop(i)};
    out = equip.outlet{on_loop(i)};
    switch equip.type{on_loop(i)}
        case 'Pipe:Adiabatic'
            node.demand_temperature(out) = node.demand_temperature(in);
        case 'Mixer'
            if node.demand_flow(out)>0
                node.demand_temperature(out) = sum(node.demand_temperature(in).*node.demand_flow(in))/node.demand_flow(out);
            else%avoid divide by zero
                node.demand_temperature(out) = mean(node.demand_temperature(in));
            end
        case 'Splitter'
            node.demand_temperature(out) = node.demand_temperature(in);
        case {'Coil:Cooling:Water';'Coil:Cooling:Water:DetailedGeometry';'Coil:Heating:Water';'Chiller:Electric:ReformulatedEIR';'Chiller:Electric:EIR';}
            if node.demand_flow(out)>0
                node.demand_temperature(out) = node.demand_temperature(in) + node.load(out)/(Cp_water*node.demand_flow(out));
            else%avoid divide by zero
                node.demand_temperature(out) = node.demand_temperature(in);
            end
        case 'WaterUse:Connections'
            we = f_index(1,strcmp(equip.name{on_loop(i)},building.water_use.name) & strcmp('Equipment',building.water_use.type));
            if isempty(building.water_use.cold_supply_temperature_schedule{we})
                node.demand_temperature(out) = T_mains;
            else
                node.demand_temperature(out) = profile.(building.water_use.cold_supply_temperature_schedule{we});
            end
        otherwise
            disp(strcat('need additional plant demand loop component of type__',equip.type{on_loop(i)}))
    end
end
end%Ends function demand_loop_energy_balance

function [plant_nodes,e_use,T_tank,tank_on] = supply_loop(building,profile,plant_nodes,e_use,T_tank,tank_on,T_air,w_air,l,dt)
Cp_water = 4186; %J/kg*K
fresh_air.T = T_air;
fresh_air.w = w_air;
s_equip = building.plant_supply_equip;
d_equip = building.plant_demand_equip;
on_loop = f_index(l,s_equip.loop);
for i = 1:1:length(on_loop)
    T_target = building.plant_loop.exit_temperature(l);
    node = building.plant_loop.setpoint_node(l);
    T_target = manager_override(building,profile,node,T_target,[],fresh_air,[]);
    in = s_equip.inlet{on_loop(i)};
    out = s_equip.outlet{on_loop(i)};
    switch s_equip.type{on_loop(i)}
        case 'Pipe:Adiabatic'
            plant_nodes.supply_temperature(out) = plant_nodes.supply_temperature(in);
        case 'Mixer'
            if plant_nodes.supply_flow(out)>0
                plant_nodes.supply_temperature(out) = sum(plant_nodes.supply_temperature(in).*plant_nodes.supply_flow(in))/plant_nodes.supply_flow(out);
            else%avoid divide by zero
                plant_nodes.supply_temperature(out) = mean(plant_nodes.supply_temperature(in));
            end
        case 'Splitter'
            plant_nodes.supply_temperature(out) = plant_nodes.supply_temperature(in);
        case {'Chiller:Electric:ReformulatedEIR';'Chiller:Electric:EIR'}
            [PLR,Q_evap,plant_nodes.supply_temperature(out),Q_avail,c] ...
                = chiller_PLR(building,s_equip.name{on_loop(i)},plant_nodes.supply_flow(in),plant_nodes.supply_temperature(in),profile);
            if PLR > 0
                cycling_ratio = min(PLR/building.chiller.min_part_load(c),1);
                PLR = max(PLR,building.chiller.min_unload_ratio(c));
                Q_false_load = Q_avail*PLR*cycling_ratio - Q_evap;
                if isfield(building.chiller,'electric_in_out_curve_type')
                    if strcmp(building.chiller.electric_in_out_curve_type{c},'LeavingCondenserWaterTemperature')
                        C_EIR_PLR = eval_curve(building.curves,building.chiller.electric_in_out_part_load_curve{c},[building.chiller.condenser_water_temperature_ref(c),PLR]);
                    else
                        dT_ref = building.chiller.condenser_water_temperature_ref(c) - building.chiller.chilled_water_temperature_ref(c);
                        dT = (building.chiller.condenser_water_temperature_ref(c) - plant_nodes.supply_temperature(out))/dT_ref;
                        T_dev = (plant_nodes.supply_temperature(out) - building.chiller.chilled_water_temperature_ref(c))/dT_ref;
                        C_EIR_PLR = eval_curve(building.curves,building.chiller.electric_in_out_part_load_curve{c},[dT,PLR,T_dev]);
                    end
                else
                    C_EIR_PLR = eval_curve(building.curves,building.chiller.electric_in_out_part_load_curve{c},PLR);
                end
                CEIR_T = eval_curve(building.curves,building.chiller.electric_in_out_temperature_curve{c},[plant_nodes.supply_temperature(out),building.chiller.condenser_water_temperature_ref(c)]);
                P_chiller = PLR*Q_avail*(1/building.chiller.COP_ref(c))*CEIR_T*C_EIR_PLR*cycling_ratio;
                Q_cond = (P_chiller*building.chiller.frac_power_2_condenser(c))+Q_evap+Q_false_load;
                if isfield(building.chiller,'condenser_type') && strcmpi(building.chiller.condenser_type(c),'aircooled')
                    %% Do something for air_cooled chiller?
                else%water cooled
                    d_e_num = f_index(s_equip.name{on_loop(i)},d_equip.name);
                    plant_nodes.load(d_equip.outlet{d_e_num}) = Q_cond;
                    plant_nodes.demand_flow(d_equip.outlet{d_e_num}) = Q_cond/(Cp_water*building.plant_loop.temperature_difference(l));
                end
                e_use.cool_elec = e_use.cool_elec + P_chiller;%Watts
            end
        case 'Boiler:HotWater'
            [PLR,load,plant_nodes.supply_temperature(out),b] = boiler_PLR(building,s_equip.name{on_loop(i)},plant_nodes.supply_flow(in),plant_nodes.supply_temperature(in),profile);
            if load>0 
                if iscell(building.boiler.normalized_eff_curve(b)) && ~isempty(building.boiler.normalized_eff_curve{b})
                    eff = eval_curve(building.curves,building.boiler.normalized_eff_curve{b},[PLR,plant_nodes.supply_temperature(out)]);
                else
                    eff = building.boiler.efficiency(b);
                end
                e_use.heat_gas = e_use.heat_gas + load/eff;%Watts
            end
        case {'Pump:VariableSpeed','Pump:ConstantSpeed'}
            %% pump heat goes into tank at start of half-loop
            plant_nodes.supply_temperature(out) = plant_nodes.supply_temperature(in);
        case 'WaterHeater:Mixed'
            wh = f_index(s_equip.name{on_loop(i)},building.water_heater.name);
            T_exit = profile.(building.water_heater.temperature_schedule{wh});
            [T_tank(wh),tank_on(wh),e_use] = water_heater(building,e_use,T_tank(wh),tank_on(wh),T_exit,T_air,plant_nodes.supply_temperature(in),plant_nodes.supply_flow(in),dt,wh);
            plant_nodes.supply_temperature(out) = T_tank(wh);
        case {'CoolingTower:SingleSpeed';'CoolingTower:TwoSpeed';'CoolingTower:VariableSpeed:Merkel'}
            [P_fan,plant_nodes.supply_temperature(out),bypass] = cooling_tower(building,s_equip.name{on_loop(i)},plant_nodes.supply_flow(in),plant_nodes.supply_temperature(in),T_target,T_air,w_air);
            e_use.tower_elec = e_use.tower_elec + P_fan;%Watts
        otherwise
            disp(strcat('need additional plant supply loop component__',s_equip.type{on_loop(i)}))
    end
end
end%Ends function supply_loop

function flow_out = load_distribution(building,out,flow_in,T_in,net_load,method,profile)
%% need equipment on the loop to add up to the load
%% equipment + bypass flows = flow in
%not doing anything if equipment is in series on a supply branch (not sure if this is ever a problem)
names = {};
bypass = [];
for j = 1:1:length(building.plant_supply_equip.name)
    if any(out == building.plant_supply_equip.inlet{j}(1))
        names(end+1,1) = building.plant_supply_equip.name(j);
        if building.plant_supply_equip.bypass(j)
            bypass = length(names);
        end
    end
end
n = length(names);
flow_out = zeros(n,1);
%first distribution of flow
switch method
    case {'OPTIMAL';'Optimal';'optimal';}
        flow_out = zeros(n,1);
        if isempty(bypass) || bypass ~=1
            flow_out(1) = flow_in;
        else
            flow_out(2) = flow_in;
        end
    case {'UNIFORMLOAD';'UniformLoad';'uniformload';}
        %split evenly amongst non-bypass branches
        if isempty(bypass)%no bypass branch
            flow_out = ones(n,1)*flow_in/n;
        else
            flow_out = ones(n,1)*flow_in/(n-1);
            flow_out(bypass) = 0;
        end
    case {'SEQUENTIALLOAD';'SequentialLoad';'sequentialload';}
        flow_out = zeros(n,1);
        if isempty(bypass) || bypass ~=1
            k = 1;
        else
            k = 2;
        end
        flow_out(k) = flow_in;
    case {'UNIFORMPLR';'UniformPLR';'uniformplr';}
        disp('need work in load_distribution of UniformPLR')
%         while abs(error_load)>tol && any(abs(error_PLR)>tol)
%         end
    case {'SEQUENTIALUNIFORMPLR';'SequentialUniformPLR';}
        disp('need work in load_distribution of SequentialUniformPLR')
%         while abs(error_load)>tol && any(abs(error_PLR)>tol)
%         end
    otherwise
        disp('need work in load_distribution of manage_plant_equip')
end

%% Re-distribution of flow
error_load = 1;
tol = 1e-2;
if net_load>0
    while abs(error_load)>tol
        [PLR,load,flow_max,PLR_opt] = PLR_load(building,names,flow_out,T_in,profile);
        flow_out = min(flow_max,flow_out);
        if isempty(load)%cooling tower, cant calculate load yet
            error_load = 0;
        else
            error_load = (net_load - sum(load))/net_load;
            sys_with_spare_cap = f_index(1,load>0 & flow_out<flow_max);
            if error_load>0 && isempty(sys_with_spare_cap)
                error_load = 0;
            end
        end
        switch method
            case {'OPTIMAL';'Optimal';'optimal';}
                if error_load>tol && any(PLR<PLR_opt)
                   %increase load up to PLR_opt in order of equipment 
                   for j = 1:1:n
                       if PLR(j)<PLR_opt(j)
                           flow_out(j) = flow_out(j)*min(PLR_opt(j)/PLR(j),(load(j)+error_load*net_load)/load(j));
                           break;
                       end
                   end
                elseif error_load<-tol && any(PLR>PLR_opt)
                    %decrease load to PLR_opt
                    for j = 1:1:n
                       if PLR(j)>PLR_opt(j)
                           flow_out(j) = flow_out(j)*min(PLR_opt(j)/PLR(j),(load(j)+error_load*net_load)/load(j));
                           break;
                       end
                   end
                elseif abs(error_load)>tol
                    %distrubte error in load uniformly
                    sys_with_spare_cap = f_index(1,load>0 & flow_out<flow_max);
                    flow_out(sys_with_spare_cap) = flow_out(sys_with_spare_cap).*(load(sys_with_spare_cap) + (net_load - sum(load))/length(sys_with_spare_cap))./load(sys_with_spare_cap);
                end
            case {'UNIFORMLOAD';'UniformLoad';'uniformload';}
                flow_out(sys_with_spare_cap) = flow_out(sys_with_spare_cap).*(load(sys_with_spare_cap) + (net_load - sum(load))/length(sys_with_spare_cap))./load(sys_with_spare_cap);
                flow_out = flow_out*flow_in/sum(flow_out); %normalize to correct total mass flow
                error_load = 0;
            case {'SEQUENTIALLOAD';'SequentialLoad';'sequentialload';}
                if error_load>tol 
                   %increase load up to max PLR in order of equipment 
                   if flow_out(k)<flow_max(k)
                       flow_out(k) = min(flow_max(k),flow_out(k)*(load(k)+error_load*net_load)/load(k));
                   else
                       k = k+1;
                       if k == bypass
                           k = k+1;
                       end
                       if k<=n
                           flow_out(k) = 0.5*flow_max(k);
                       end
                   end
                elseif error_load<-tol
                    %decrease load 
                    flow_out(k) = max(0,flow_out(k)*(load(k)+error_load*net_load)/load(k));
                    if flow_out(k) == 0
                        k = k-1;
                        if bypass == k
                            k = k-1;
                        end
                    end
                end
            case {'UNIFORMPLR';'UniformPLR';'uniformplr';}
                disp('need work in load_distribution of UniformPLR')
        %         while abs(error_load)>tol && any(abs(error_PLR)>tol)
        %         end
            case {'SEQUENTIALUNIFORMPLR';'SequentialUniformPLR';}
                disp('need work in load_distribution of SequentialUniformPLR')
        %         while abs(error_load)>tol && any(abs(error_PLR)>tol)
        %         end
            otherwise
                disp('need work in load_distribution of manage_plant_equip')
        end
        if ~isempty(bypass)
            flow_out(bypass) = max(0,flow_in - (sum(flow_out)- flow_out(bypass)));
        end
    end
end
end%Ends function load_distribution

function [PLR,load,flow_max,PLR_opt] = PLR_load(building,names,flow,T_in,profile)
Cp_water = 4186; %J/kg*K
n = length(names);
PLR = zeros(n,1);
load = zeros(n,1);
flow_max = zeros(length(names),1);
PLR_opt = zeros(n,1);
for j = 1:1:n
    k = f_index(names{j},building.plant_supply_equip.name);
    if any(strcmp(building.plant_supply_equip.type{k},{'Chiller:Electric:ReformulatedEIR';'Chiller:Electric:EIR';'Chiller:ConstantCOP';}))
        [PLR(j),load(j),T_out,Q_avail,c] = chiller_PLR(building,names{j},flow(j),T_in,profile); 
        if PLR(j)*Q_avail/(Cp_water*(T_in - T_out))<flow(j)
            flow_max(j) = PLR(j)*Q_avail/(Cp_water*(T_in - T_out));
        else
            flow_max(j) = building.chiller.chilled_water_flow_ref(c)*building.plant_loop.fluid_density(building.plant_supply_equip.loop(k));%max flow in kg/s
        end
        PLR_opt(j) = building.chiller.opt_part_load(c);
    elseif any(strcmp(building.plant_supply_equip.type{k},{'Boiler:HotWater'}))
        [PLR(j),load(j),T_out,b] = boiler_PLR(building,names{j},flow(j),T_in,profile);
        if load(j)/(Cp_water*(T_out - T_in))<flow(j)
            flow_max(j) = load(j)/(Cp_water*(T_out - T_in));
        else
            flow_max(j) = building.boiler.design_flow_rate(b)*building.plant_loop.fluid_density(building.plant_supply_equip.loop(k));%max flow in kg/s
        end
        PLR_opt(j) = building.boiler.optimal_part_load(b);
    elseif any(strcmp(building.plant_supply_equip.type{k},{'CoolingTower:SingleSpeed'}))
        ct = f_index(names{j},building.cooling_tower.name);
        load = [];
        PLR(j) = 1;
        PLR_opt(j) = 1;
        flow_max(j) = building.cooling_tower.water_flow(ct)*building.plant_loop.fluid_density(building.plant_supply_equip.loop(k));
    elseif any(strcmp(building.plant_supply_equip.type{k},{'WaterHeater:Mixed'}))
        load = [];
        PLR(j) = 1;
        PLR_opt(j) = 1;
        flow_max(j) = building.plant_loop.max_flow_rate(building.plant_supply_equip.loop(k));
    elseif any(strcmp(building.plant_supply_equip.type{k},{'Chiller:EngineDriven';'Chiller:Absorption';'Chiller:Absorption:Indirect';'Chiller:CombustionTurbine';'ChillerHeater:Absorption:DirectFired';'ChillerHeater:Absorption:DoubleEffect';'CoolingTower:TwoSpeed';'CoolingTower:VariableSpeed:Merkel'}))
            disp(strcat('need examples of these types of equipment__',building.plant_supply_equip.type{k}))
    elseif  building.plant_supply_equip.bypass(k)
        flow_max(j) = sum(flow);
    end
end
end%Ends function PLR_load

function [PLR,load,T_out,b] = boiler_PLR(building,b_name,flow,T_in,profile)
Cp_water = 4186; %J/kg*K
b = f_index(b_name,building.boiler.name);
if strcmp(building.boiler.flow_mode{b},'LeavingSetpointModulated')
    m = f_index(building.boiler.water_outlet{b},building.manager.node);
    T_out = profile.(building.manager.schedule{m});
else 
    disp('need some other way to control boiler load')
    T_out = building.boiler.design_water_temperature(b);
end
load = flow*Cp_water*(T_out- T_in);
PLR = min(load/building.boiler.capacity(b),building.boiler.max_part_load(b));
if load>0 && PLR<building.boiler.min_part_load(b)
    PLR = building.boiler.min_part_load(b);
    load = PLR*building.boiler.capacity(b);
end
end%Ends function boiler_PLR

function [PLR,load,T_out,Q_avail,c] = chiller_PLR(building,c_name,flow,T_in,profile)
%%see pg 763 about standard IPLV values
Cp_water = 4186; %J/kg*K
c = f_index(c_name,building.chiller.name);
if strcmp(building.chiller.flow_mode{c},'LeavingSetpointModulated')
    m = f_index(building.chiller.chilled_water_outlet{c},building.manager.node);
    T_set_out = profile.(building.manager.schedule{m});
else 
    disp('need some other way to control chiller load')
    T_set_out = building.chiller.chilled_water_temperature_ref(c);
end
Q_evap = flow*Cp_water*(T_in - T_set_out);% (Watts) kg/s*J/kgK * K = J/s
%need to update with actual condenser water outlet temperature
CCFT = eval_curve(building.curves,building.chiller.capacity_temperature_curve{c},[T_set_out,building.chiller.condenser_water_temperature_ref(c)]);
Q_avail = building.chiller.capacity(c)*CCFT;%(Watts)
min_T_out = T_in - Q_avail/(Cp_water*flow);
T_out = max(T_set_out,min_T_out);
PLR = min(Q_evap/Q_avail,building.chiller.max_part_load(c));
load = PLR*Q_avail;
end%Ends function chiller_PLR