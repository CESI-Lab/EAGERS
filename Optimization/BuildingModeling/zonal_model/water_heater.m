function [T,tank_on,e_use] = water_heater(building,e_use,T,tank_on,T_exit,T_air,T_source,flow_source,dt,wh)
%%not sure what to do with use side heat exchanger of tank
Cp_water = 4186; %J/kg*K
den_water = 997;%kg/m^3

e_source = building.water_heater.source_side_effectiveness(wh);
m_cp = building.water_heater.tank_volume(wh)*den_water*Cp_water;
T_dif = building.water_heater.deadband_temperature_dif(wh);
if isnan(T_dif)
    T_dif = 0;
end
if isnan(e_source)
    e_source = 1;
end
if isnan(building.water_heater.on_cycle_fuel_use(wh))
    q_on_cycle = 0;
else
    q_on_cycle = building.water_heater.on_cycle_fuel_use(wh)*building.water_heater.on_cycle_frac2tank(wh);
end
if isnan(building.water_heater.off_cycle_fuel_use(wh))
    q_off_cycle = 0;
else
    q_off_cycle = building.water_heater.off_cycle_fuel_use(wh)*building.water_heater.off_cycle_frac2tank(wh);
end
UA_on_cycle = building.water_heater.on_cycle_loss2ambient(wh);
UA_off_cycle = building.water_heater.off_cycle_loss2ambient(wh);
q_heater = building.water_heater.heater_maximum_capacity(wh)*building.water_heater.heater_efficiency(wh);
a_on = (q_heater + q_on_cycle + UA_on_cycle*T_air + e_source*flow_source*Cp_water*T_source)/m_cp;
b_on = -(UA_on_cycle + e_source*flow_source*Cp_water)/m_cp;
a_off = (q_off_cycle + UA_off_cycle*T_air + e_source*flow_source*Cp_water*T_source)/m_cp;
b_off = -(UA_off_cycle + e_source*flow_source*Cp_water)/m_cp;

sec_on = 0;
t = 0;
while t < dt
    step = min(60,dt-t);%minutes
    if T<T_exit-T_dif
        tank_on = true;
    end
    if tank_on
        if b_on == 0
            T_new = a_on*step + T;
        else
            T_new = (a_on/b_on + T)*exp(b_on*step) - a_on/b_on;
        end
        if T_new>=T_exit%turns off partway through step
            tank_on = false;
            if b_on == 0
                t_on = (T_exit - T)/a_on;
            else
                t_on = 1/b_on*log((a_on/b_on + T_exit)/(a_on/b_on + T));
            end
            sec_on = sec_on + t_on;%on for fraction of minute
            %off for remainder of minute
            T = T_exit;
            if b_off == 0
                T = a_off*(step-sec_on) + T;
            else
                T = (a_off/b_off + T)*exp(b_off*(step-sec_on)) - a_off/b_off;
            end
        else
            sec_on = sec_on + step;
        end
    else
        if b_off == 0
            T = a_off*step + T;
        else
            T = (a_off/b_off + T)*exp(b_off*step) - a_off/b_off;
        end
    end
    t = t + step;
end

if strcmp(building.water_heater.heater_fuel(wh),'NATURALGAS')
    e_use.water_gas = e_use.water_gas + sec_on/dt*building.water_heater.heater_maximum_capacity(wh);%Watts
else
    e_use.water_elec = e_use.water_elec + sec_on/dt*building.water_heater.heater_maximum_capacity(wh);%Watts
end
if q_on_cycle>0
    if strcmp(building.water_heater.on_cycle_fuel_type(wh),'NATURALGAS')
        e_use.water_gas = e_use.water_gas + sec_on/dt*building.water_heater.on_cycle_fuel_use(wh);%Watts
    else
        e_use.water_elec = e_use.water_elec + sec_on/dt*building.water_heater.on_cycle_fuel_use(wh);%Watts
    end
end
if q_off_cycle>0
    if strcmp(building.water_heater.off_cycle_fuel_type(wh),'NATURALGAS')
        e_use.water_gas = e_use.water_gas + (dt-sec_on)/dt*building.water_heater.off_cycle_fuel_use(wh);%Watts
    else
        e_use.water_elec = e_use.water_elec + (dt-sec_on)/dt*building.water_heater.off_cycle_fuel_use(wh);%Watts
    end
end
end%Ends function water_heater