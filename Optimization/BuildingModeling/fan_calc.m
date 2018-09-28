function varargout = fan_calc(varargin)
fans = varargin{1};
name = varargin{2};
k = nonzeros((1:length(fans.name))'.*strcmpi(name,fans.name));
mode = varargin{3};
air_density = 1.225; %kg/m^3
switch mode
    case 'sizing'
        %do nothing
        
    case 'operation'
        air = varargin{4};
        profile = varargin{5};
        flow = min(air.m_dot,fans.flow_rate(k)*profile);
        Q_tot = flow.*fans.pressure_rise(k)/(fans.fan_efficiency(k)*air_density);% kg/s * Pa / (kg/m^3) = N*m/s = W
        Q_shaft = fans.motor_efficiency(k)*Q_tot;
        Q_to_air = Q_shaft + (Q_tot-Q_shaft)*fans.motor_frac(k);
        air_out.m_dot = flow;
        air_out.w = air.w;
        if air.m_dot>0
            air_out.h = air.h + Q_to_air/1000/air.m_dot;%all h in kJ/kg
        else
            air_out.h = air.h;
        end
        air_out.T = (air_out.h - air_out.w*2501)./(1.006 + air_out.w*1.86);%Temperature in C using enthalpy in kJ/kg and humidity in massH2O/massAir
        varargout = {air_out,Q_tot,Q_to_air};
end