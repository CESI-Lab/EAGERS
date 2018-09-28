function varargout = heating_coil(varargin)
coil = varargin{1};
name = varargin{2};
k = nonzeros((1:length(coil.name))'.*strcmpi(name,coil.name));
mode = varargin{3};
switch mode
    case 'sizing'
        switch coil.type{k}
            case 'Water'
                Tw_in0 = coil.inlet_water_temperature(k);
                air_in.T = varargin{4};
                air_out.T = coil.air_outlet_temperature(k);
                air_in.w = varargin{5};
                air_out.w = varargin{5};%%assume same humidity
                air_in.h = psychometric_h(air_in);
                air_out.h = psychometric_h(air_out);
                m_a = varargin{6};
                Cp_air = 1006 + 1860*air_in.w; %J/kg*K
                Cp_water = 4186; %J/kg*K
                coil.capacity(k,1) = m_a*(air_out.h-air_in.h)*1000;%W
                m_w = coil.capacity(k,1)/(Cp_water*(coil.inlet_water_temperature(k) - coil.outlet_water_temperature(k)));
                
                C_air = Cp_air*m_a;
                C_w = Cp_water*m_w;
                C_min = min(C_air,C_w);
                Z = C_min/max(C_air,C_w);
                max_heat_transfer = C_min*(Tw_in0 - air_in.T);
                effectiveness = (air_out.T - air_in.T)*C_air/max_heat_transfer;
                NTU = [0,.001,.01,.1,1,logspace(1,5,5)];
                effect = 1 - exp((exp(-Z*NTU.^.78)-1)./(Z*NTU.^(-.22)));
                for j = 1:1:3
                    nn = nnz(effect<effectiveness);
                    NTU = linspace(NTU(nn),NTU(11-nn),10); 
                    effect = 1 - exp((exp(-Z*NTU.^.78)-1)./(Z*NTU.^(-.22)));
                end
                NTU = interp1(effect,NTU,effectiveness);
                coil.UA(k,1) = NTU*C_min;
                coil.rated_air_flow(k,1) = m_a;
                coil.rated_water_flow(k,1) = m_w;
                if isnan(coil.air_inlet_temperature(k,1))
                    coil.air_inlet_temperature(k,1) = air_in.T;
                end
                varargout = {coil};
            case {'Fuel';'Electric'}
                air_in = varargin{4};
                T_supply = varargin{5};
                m_a = varargin{6};
                coil.capacity(k,1) = varargin{7};
                varargout = {coil};
            otherwise
                disp('need to add missing heating coil type')
        end
    case 'operation'
        air_in = varargin{4};
        sensible_load = varargin{5};
        air_out = air_in;
        switch coil.type{k}
            case {'Fuel';'Electric'}
                if air_in.m_dot>1e-8
                    air_out.h = (air_in.h.*air_in.m_dot + sensible_load/1000)/air_out.m_dot;
                    air_out.w = air_in.w; %humidity
                    air_out.T = (air_out.h - air_out.w*2501)./(1.006 + air_out.w*1.86);%Temperature in C using enthalpy in kJ/kg and humidity in massH2O/massAir
                end
                heat = sensible_load/coil.efficiency(k); %heat added in W
                varargout = {air_out,heat};
            case 'Water'
                Tw_in = varargin{6}; %supply water
                Cp_water = 4186; %J/kg*K
                m_w = max(0,sensible_load)/(Cp_water*(Tw_in - coil.outlet_water_temperature(k)));
                Cp_air = 1006 + 1860*air_in.w; %J/kg*K
                Tw_in0 = coil.inlet_water_temperature(k);
                Ta_in0 = coil.air_inlet_temperature(k);
                Ta_in = air_in.T;
                m_a = air_in.m_dot;
                
                if sensible_load>0 %heating
                    %part of UA calculations depending only on inputs
                    if ~isempty(coil.air_water_convection_ratio{k})
                        r = str2double(coil.air_water_convection_ratio{k});
                    else
                        r = 1;
                    end
                    hA_w0 = coil.UA(k)*(r+1/r);
                    nfhA_a0 = r*hA_w0;
                    x_a = 1 + 4.769e-3*(Ta_in - Ta_in0);
                    nfhA_a = x_a*((m_a/coil.rated_air_flow(k))^.8)*nfhA_a0;
                    x_w = 1 + (.014/(1 + .014*Tw_in0))*(Tw_in - Tw_in0);
                    
                    %%loop to solve for water flow that gives correct sensible load (page 960 of the manual)
                    %%water flow rate and water outlet temperature change
                    tolerance = 1e-4;
                    error = 1;
                    while abs(error)>tolerance
                        %find new UA factor
                        hA_w = x_w*((m_w/coil.rated_water_flow(k))^.85)*hA_w0;
                        UA = 1/(1/hA_w + 1/nfhA_a);

                        %capacitance
                        C_air = Cp_air*m_a;
                        C_w = Cp_water*m_w;
                        C_min = min(C_air,C_w);
                        Z = C_min/max(C_air,C_w);

                        %Find outlet temperatures
                        NTU = UA/C_min;
                        effect = 1 - exp((exp(-Z*NTU.^.78)-1)./(Z*NTU.^(-.22)));
                        Ta_out = Ta_in + effect*C_min*(Tw_in - Ta_in)/C_air;
                        
                        air_out.T = Ta_out;
                        air_out.h = psychometric_h(air_out);
                        actual_load = air_out.m_dot*(air_out.h - air_in.h)*1000;
                        error = (sensible_load - actual_load)/max(sensible_load,actual_load);
                        m_w = min(coil.rated_water_flow(k),m_w*(1+error));
                        if error>0 && coil.rated_water_flow(k)==m_w%exits the loop when it reaches max water flow
                            error = 0;
                        end
                    end
                    water_out.m_dot = m_w;
                    water_out.T = Tw_in - C_air*(Ta_out - Ta_in)/C_w;
%                     heat = C_w*(Tw_in - water_out.T);%Q_coil (W)
                else
                    water_out.m_dot = 0;
                    water_out.T = coil.outlet_water_temperature(k);
                end
                varargout = {air_out,water_out};
            otherwise
                disp('need to add missing heating coil type')
        end
end
end%Ends function heating_coil