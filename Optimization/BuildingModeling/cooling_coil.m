function varargout = cooling_coil(varargin)
coil = varargin{1};
k = nonzeros((1:length(coil.name))'.*strcmpi(varargin{2},coil.name));
mode = varargin{3};
switch mode
    case 'sizing'
        switch coil.type{k}
            case {'DX:TwoSpeed';'DX:SingleSpeed';}%page 805
                if isnan(coil.rated_capacity(k)) %needs sizing
                    curve = varargin{4};
                    flow = varargin{5};
                    T_sup = varargin{6};
                    w_sup = varargin{7};
                    all_outdoor = varargin{8};
                    coil.rated_air_flow(k,1) = flow;
                    T_amb = 35;
                    if all_outdoor
                        T_mix = 26.7;
                    else
                        T_mix = 35;
                    end
                    w_mix = 0.01115; %T_wb = 19.5, RH = 50.7%
                    h_mix = (1.006*T_mix + w_mix.*(2501 + 1.86*T_mix));%kJ/kg %ASHRAE 2013 fundamentals eq. 32 (http://www.ce.utexas.edu/prof/Novoselac/classes/ARE383/Handouts/F01_06SI.pdf)

                    h_sup = (1.006*T_sup + w_sup.*(2501 + 1.86*T_sup));%kJ/kg %ASHRAE 2013 fundamentals eq. 32 (http://www.ce.utexas.edu/prof/Novoselac/classes/ARE383/Handouts/F01_06SI.pdf)
                    peak_capacity = flow*(h_mix-h_sup)*1000;% kg/s * kJ/kg * 1000 J/kJ = W
                    CapMod = eval_curve(curve,coil.capacity_v_temperature_curve{k},[T_mix,T_amb]); 
                    CC_rated = peak_capacity/CapMod;

                    if all_outdoor
                        CC_rated = max(min(CC_rated,flow/0.00001677),flow/0.00003355);
                    else
                        CC_rated = max(min(CC_rated,flow/0.00004027),flow/0.00006041);
                    end
                    coil.rated_capacity(k,1) = CC_rated;
                    switch coil.type{k}
                        case 'DX:TwoSpeed'
                            coil.rated_capacity2(k,1) = CC_rated/2;
                            coil.rated_air_flow2(k,1) = flow/2;
                            varargout = {coil};
                        case 'DX:SingleSpeed'
                            %%%
                    end
                end
                varargout = {coil};
            case 'Water:Detailed'
                air_in.T = varargin{4};
                air_out.T = coil.air_outlet_temperature(k);
                air_in.w = varargin{5};
                air_out.w = varargin{5};%%assume same humidity
                air_in.h = psychometric_h(air_in);
                air_out.h = psychometric_h(air_out);
                m_a = varargin{6};
                Cp_water = 4186; %J/kg*K
                coil.capacity(k,1) = m_a*(air_in.h-air_out.h)*1000;%W
                m_w = coil.capacity(k,1)/(Cp_water*(coil.outlet_water_temperature(k) - coil.inlet_water_temperature(k)));
                coil.tubes_per_row(k,1) = max(3,13750*m_w);
                coil.fin_diameter(k,1) = 0.335*m_a;
                coil.airflow_area(k,1) = 0.44*m_a;
                coil.fin_area(k,1) = 78.5*m_a;
                coil.area_inside(k,1) = 4.4*coil.tube_id(k)*coil.tube_rows(k)*coil.tubes_per_row(k,1);
                coil.area_outside(k,1) = 4.1*coil.tube_od(k)*coil.tube_rows(k)*coil.tubes_per_row(k,1);
                coil.coil_depth(k,1) = coil.tube_spacing(k)*coil.tube_rows(k);
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
                coil.capacity(k,1) = m_a*(air_in.h-air_out.h)*1000;%W
                m_w = coil.capacity(k,1)/(Cp_water*(coil.outlet_water_temperature(k) - coil.inlet_water_temperature(k)));
                %find nominal UA factor
                C_air = Cp_air*m_a;
                C_w = Cp_water*m_w;
                C_min = min(C_air,C_w);
                Z = C_min/max(C_air,C_w);
                max_heat_transfer = C_min*(air_in.T - Tw_in0);
                effectiveness = (air_in.T - air_out.T)*C_air/max_heat_transfer;
                NTU = [0,.001,.01,.1,1,logspace(1,5,5)];
                effect = 1 - exp((exp(-Z*NTU.^.78)-1)./(Z*NTU.^(-.22)));
                for j = 1:1:3
                    NTU = linspace(NTU(nnz(effect<effectiveness)),NTU(11-nnz(effect>effectiveness)),10); 
                    effect = 1 - exp((exp(-Z*NTU.^.78)-1)./(Z*NTU.^(-.22)));
                end
                NTU = interp1(effect,NTU,effectiveness);
                coil.UA(k,1) = NTU*C_min;
                coil.rated_air_flow(k,1) = m_a;
                coil.rated_water_flow(k,1) = m_w;
                coil.air_inlet_temperature(k,1) = air_in.T;
                varargout = {coil};
            otherwise
                disp('need to add missing cooling coil type')
        end
    case 'operation'
        switch coil.type{k}
            case {'DX:TwoSpeed';'DX:SingleSpeed';}%page 805
                curve = varargin{4};
                air_in = varargin{5};
                T_db = varargin{6}; %outside dry air temperature
                sensible_load = varargin{7};

                air_out = air_in;
                if air_out.m_dot>1e-8 && sensible_load>0
                    air_out.h = (air_in.h.*air_in.m_dot - sensible_load/1000)/air_out.m_dot;
                    air_out.w = ones(length(air_in.T),1)*varargin{8}; %humidity
                    air_out.T = (air_out.h - air_out.w*2501)./(1.006 + air_out.w*1.86);%Temperature in C using enthalpy in kJ/kg and humidity in massH2O/massAir

                    %add details if exceeding max load, or if water condenses
                    switch coil.type{k}
                        case 'DX:TwoSpeed'
                            input(:,1) = air_in.T; %needs to be wet-bulb, once I figure out humidity
                            input(:,2) = T_db; %needs to be outdoor wetbulb temperature if it is an evaporative cooled condenser
                            TotCapTempMod_1 = eval_curve(curve,coil.capacity_v_temperature_curve{k},input);
                            EIRTempMod_1 = eval_curve(curve,coil.energy_input_v_temperature_curve{k},input);
                            TotCapTempMod_2 = eval_curve(curve,coil.capacity_v_temperature_curve2{k},input);
                            EIRTempMod_2 = eval_curve(curve,coil.energy_input_v_temperature_curve2{k},input);

                            input2(:,1) = air_in.m_dot/coil.rated_air_flow(k); %normalized flow
                            TotCapFlowMod = eval_curve(curve,coil.capacity_v_flow_curve{k},input2); 
                            EIRFlowMod = eval_curve(curve,coil.energy_input_v_flow_curve{k},input2);  

                            max_load_1 = coil.rated_capacity(k)*TotCapTempMod_1*TotCapFlowMod;
                            max_load_2 = coil.rated_capacity2(k)*TotCapTempMod_2*TotCapFlowMod;
                            EIR_1 = (1/coil.rated_COP(k))*EIRTempMod_1*EIRFlowMod;%energy input ratio
                            EIR_2 = (1/coil.rated_COP2(k))*EIRTempMod_2*EIRFlowMod;%energy input ratio

                            Speed_ratio = max(0,(sensible_load - max_load_2)/(max_load_1 - max_load_2));
                            if Speed_ratio>0
                                Power = Speed_ratio*max_load_1*EIR_1 + (1-Speed_ratio)*max_load_2*EIR_2;
                                heat_rejected = Speed_ratio*max_load_1*(1+EIR_1) + (1-Speed_ratio)*max_load_2*(1+EIR_2);
                            else
                                PLR = sensible_load/(coil.rated_capacity(k)*coil.rated_sensible_heat_ratio(k)); %
                                PartLoadFrac = max(.7,eval_curve(curve,coil.part_load_curve{k},PLR));
                                RTF = PLR./PartLoadFrac;%runtime fraction
                                Power = max_load_2*EIR_2*RTF;
                                heat_rejected = max_load_2*(1+EIR_2);
                            end  
                        case 'DX:SingleSpeed'
                            input(:,1) = air_in.T; %needs to be wet-bulb, once I figure out humidity
                            input(:,2) = T_db; %needs to be outdoor wetbulb temperature if it is an evaporative cooled condenser
                            TotCapTempMod= eval_curve(curve,coil.capacity_v_temperature_curve{k},input);
                            EIRTempMod = eval_curve(curve,coil.energy_input_v_temperature_curve{k},input);

                            input2(:,1) = air_in.m_dot/coil.rated_air_flow(k); %normalized flow
                            TotCapFlowMod = eval_curve(curve,coil.capacity_v_flow_curve{k},input2); 
                            EIRFlowMod= eval_curve(curve,coil.energy_input_v_flow_curve{k},input2);  
                            PLR = sensible_load/(coil.rated_capacity(k)*coil.rated_sensible_heat_ratio(k)); %
                            PartLoadFrac = max(.7,eval_curve(curve,coil.part_load_curve{k},PLR));

                            max_load = coil.rated_capacity(k)*TotCapTempMod*TotCapFlowMod;
                            EIR = (1/coil.rated_COP(k))*EIRTempMod*EIRFlowMod;%energy input ratio
                            RTF = PLR./PartLoadFrac;%runtime fraction
                            Power = max_load*EIR*RTF;
                            heat_rejected = max_load*(1+EIR);
                    end
                    %%need to add sensible/latent heat calculations for moisture.
                else
                    Power = 0;
                end
                varargout = {air_out,Power};
            case 'Water:Detailed'
                %% see page 796 of manual. 3 air temperature states and 3 water temp states. Move boundary(middle state) until air temperature = dewpoint
                air_in = varargin{4};
                sensible_load = varargin{5};
                Tw_in = varargin{6}; %supply water
                Cp_water = 4186; %J/kg*K
                air_out = air_in;
                cool.water = 0;
                cool.electricity = 0;
                if sensible_load>0 %cooling
                    Pr = 0.733;
                    viscosity = 1.846e-5;
                    tube_thickness = coil.tube_od(k)-coil.tube_id(k);
                    %%compute dew point and humidity ratio                    
                    P = 101.325; % atmospheric pressure (kPa)
                    P_H2O = air_in.w/0.621945*P/(1+air_in.w/0.621945);
%                     Tdb_K = air_in.T+273.15; %Tdb (Kelvin)
%                     satP = exp((-5.8002206e3)./Tdb_K + 1.3914993 - 4.8640239e-2*Tdb_K + 4.1764768e-5*Tdb_K.^2 - 1.4452093e-8*Tdb_K.^3 + 6.5459673*log(Tdb_K))/1000; %saturated water vapor pressure ASHRAE 2013 fundamentals eq. 6 in kPa valid for 0 to 200C
%                     rel_humid = P_h2O/satP;
                    T_dp = 6.54 + 14.526*log(P_H2O) + 0.7389*log(P_H2O).^2 + 0.09486*log(P_H2O).^3 + 0.4569*(P_H2O).^0.1984; %Dew point from partial pressure of water using ASHRAE 2013 Fundamentals eqn 39 valid from 0C to 93C

                    Re = 4*tube_thickness*(1+air_in.w)*air_in.m_dot/(A_s_total*viscosity);
                    
                    m_w = max(0,sensible_load)/(Cp_water*(coil.outlet_water_temperature(k) - Tw_in));
                    
                    
                    
                    water_out.m_dot = m_w;
                    water_out.T = Tw_in - C_air*(Ta_out - Ta_in)/C_w;
                    cool.water = C_w*(water_out.T - Tw_in);%Q_coil (W)
                else
                    water_out.m_dot = 0;
                    water_out.T = coil.outlet_water_temperature(k);
                end
                varargout = {air_out,water_out,cool};
            case 'Water'
                %The option to identify which mode of operation the Simple mode analysis should perform ie, for
                %a given set of inputs would the coil be Dry or Wet is decided by set of conditions described below.
                %• IF (Temperature Dewpoint Air < Water Inlet Temperature) THEN the coil is Dry
                %and we call the Subroutine Coil Completely Dry. In this case outlet temperature of air would
                %be higher than the air dewpoint and hence there would be no condensation.
                %• IF (Temperature Dewpoint Air > Water Inlet Temperature) THEN the coil is
                %completely wet, call subroutine Coil Completely Wet, it is assumed that moisture condensation
                %occurs over completely surface of the coil. 
                
                air_in = varargin{4};
                air_out = air_in;
                sensible_load = varargin{5};
                Tw_in = varargin{6}; %supply water
                Cp_water = 4186; %J/kg*K
                m_w = max(0,sensible_load)/(Cp_water*(coil.outlet_water_temperature(k) - Tw_in));
                Cp_air = 1006 + 1860*air_in.w; %J/kg*K          
                Tw_in0 = coil.inlet_water_temperature(k);
                Ta_in0 = coil.air_inlet_temperature(k);
                Ta_in = air_in.T;
                m_a = air_in.m_dot;
                cool.water = 0;
                cool.electricity = 0;
                
                if sensible_load>0 %cooling
                    %part of UA calculations depending only on inputs
                    hA_w0 = coil.UA(k);
                    nfhA_a0 = hA_w0;
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
                        actual_load = air_out.m_dot*(air_in.h - air_out.h)*1000;
                        error = (sensible_load - actual_load)/max(sensible_load,actual_load);
                        m_w = min(coil.rated_water_flow(k),m_w*(1+error));
                        if error>0 && coil.rated_water_flow(k)==m_w%exits the loop when it reaches max water flow
                            error = 0;
                        end
                    end
                    water_out.m_dot = m_w;
                    water_out.T = Tw_in - C_air*(Ta_out - Ta_in)/C_w;
                    cool.water = C_w*(water_out.T - Tw_in);%Q_coil (W)
                else
                    water_out.m_dot = 0;
                    water_out.T = coil.outlet_water_temperature(k);
                end
                
                varargout = {air_out,water_out,cool};
            otherwise
                disp('need to add missing cooling coil type')
        end
end
end%Ends function cooling_coil

function out = eval_curve(curve,name,input)
k = nonzeros((1:length(curve.name))'.*strcmpi(name,curve.name));
x = min(curve.max_x(k),max(curve.min_x(k),input(:,1)));
switch curve.type{k}
    case 'Quadratic'
        out = curve.a0(k) + curve.a1(k)*x + curve.a2(k)*x.^2;
    case 'Cubic'
        out = curve.a0(k) + curve.a1(k)*x + curve.a2(k)*x.^2 + curve.a3(k)*x.^3;
    case 'Biquadratic'
        y = min(curve.max_y(k),max(curve.min_y(k),input(:,2)));
        out = curve.a0(k) + curve.a1(k)*x + curve.a2(k)*x.^2 + curve.b1(k)*y + curve.b2(k)*y.^2 + curve.ab(k)*x.*y;
    case 'Bicubic'
        y = min(curve.max_y(k),max(curve.min_y(k),input(:,2)));
        out = curve.a0(k) + curve.a1(k)*x + curve.a2(k)*x.^2 + curve.a3(k)*x.^3 + curve.b1(k)*y + curve.b2(k)*y.^2 + curve.b3(k)*y.^3 + curve.ab(k)*x.*y + curve.aab(k)*x.^2.*y + curve.abb(k)*x.*y.^2;
end
end%Ends function eval_curve