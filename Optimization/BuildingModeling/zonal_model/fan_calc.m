function [air_out,Q_tot,Q_to_air] = fan_calc(fans,name,air,avail)
k = nonzeros((1:length(fans.name))'.*strcmpi(name,fans.name));
air_density = 1.204; %kg/m^3
flow = min(air.m_dot/air_density,fans.flow_rate(k)*avail);
if flow>0
    switch fans.type{k}
        case {'ConstantVolume';'OnOff'}
            PLF = 1;%should be determined by partial loading of other components in HVAC loop, see page 1040
            RTF = flow/fans.flow_rate(k)/PLF;
            Q_tot = RTF*fans.flow_rate(k)*fans.pressure_rise(k)/fans.fan_efficiency(k);% m^3/s * Pa = N*m/s = W
        case 'VariableVolume'
            f_flow = flow/fans.flow_rate(k);
            f_pl = fans.a1(k) + fans.a2(k)*f_flow + fans.a3(k)*f_flow^2 + fans.a4(k)*f_flow^3 + fans.a5(k)*f_flow^4;
            Q_tot = f_pl*fans.flow_rate(k)*fans.pressure_rise(k)/fans.fan_efficiency(k);% m^3/s * Pa = N*m/s = W
        case 'SystemModel' %designed to replace other types
            
    end
    Q_shaft = fans.motor_efficiency(k)*Q_tot;
    Q_to_air = Q_shaft + (Q_tot-Q_shaft)*fans.motor_frac(k);%W
    cp = (1006 + 1860*air.w)*air_density;
    air_out.T = air.T + Q_to_air/(cp*flow);
    air_out.w = air.w;
    air_out.h = air.h + Q_to_air/(flow*air_density);%all h in J/kg
    air_out.flow = flow;
    air_out.m_dot = flow*air_density;
else
    air_out = air;
    Q_tot = 0;
    Q_to_air = 0;
end
end%Ends function fan_calc