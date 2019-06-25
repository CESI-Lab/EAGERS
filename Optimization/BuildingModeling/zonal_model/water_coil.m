function [air_out,m_w,actual_load,out] = water_coil(building,name,air_in,set_point,plant_nodes)
%% mode 1 (empty plant_nodes, only heating coils): Given inlet temperatures and zone load request, for given air_flow find water flow that meets load. If water flow saturates, increase air flow until load is met or air flow saturates
%% mode 2: Given water and air inlet flows and temperature, find air outlet
%% Detailed water coil applies only to cooling, only mode 2
k = f_index(name,building.coils.heating.name);
if isempty(k)
    k = f_index(name,building.coils.cooling.name);
    mode = 'cooling';
    coil = building.coils.cooling;
else
    mode = 'heating';
    coil  = building.coils.heating;
end
Cp_water = 4186; %J/kg*K
if isempty(plant_nodes) 
    sensible_load = min(set_point,building.coils.heating.capacity(k));
    air.T  = coil.air_inlet_temperature(k);
    air.w  = coil.air_inlet_humidity_ratio(k);
    air.m_dot = air_in;
    air.h = psychometric(air,'h');
    air_in = air;
    Cp_air = 1006 + 1860*air_in.w; %J/kg*K
    Tw_in = building.coils.heating.water_inlet_temperature(k);
else
    T_supply = set_point;
    n = f_index(name,building.plant_demand_equip.name);
    in = building.plant_demand_equip.inlet{n};
    out = building.plant_demand_equip.outlet{n};
    Tw_in = plant_nodes.demand_temperature(in);
    m_w = plant_nodes.demand_flow(in);
    Cp_air = 1006 + 1860*air_in.w; %J/kg*K
    switch mode
        case 'heating'
            sensible_load  = air_in.m_dot*Cp_air*(T_supply - air_in.T);
        case 'cooling'
            sensible_load  = air_in.m_dot*Cp_air*(air_in.T - T_supply);
    end
end
air_out = air_in;
if sensible_load/air_in.m_dot>1    
    tolerance = 1e-4;
    switch coil.type{k}
        case 'Water:DetailedGeometry'
            %% see page 796 of manual. 3 air temperature states and 3 water temp states. Move boundary(middle state) until air temperature = dewpoint
            %% Major problems: 1) assuming constant fin efficiency of 0.75 (cant get Bessel functions to give reasonable answer)
            %%                 2) dividing fin diameter by 100. Otherwise it was nit giving reasonable answers for resistance terms
            %%                 3) removed fouling factor for R_mf resistance. Was too high and overly limited heat transfer
            den_water = 995;%kg/m^3
            m_w_rated = coil.max_flow(k)*den_water;%convert m^3/s to kg/s
            h1 = psychometric('Twb',Tw_in,'h');
            h2 = psychometric('Twb',Tw_in+1,'h');
            bb = h2-h1;
            aa = h1 - bb*Tw_in;
            Pr = 0.733;
            viscosity = 1.846e-5;
            tube_thickness = (coil.tube_od(k)-coil.tube_id(k))/2;
            pipe_area = pi()*(coil.tube_id(k)/2)^2;%*coil.tube_rows(k)*coil.tubes_per_row(k)
            surface_area = coil.area_outside(k) + coil.fin_area(k);
            coil.fin_diameter(k) = coil.fin_diameter(k)/100;%% seems to be necessary, fin diameter is 6-16 m?
            fin_height = (coil.fin_diameter(k)-coil.tube_od(k))/2;
            
            Re_a = 4*coil.coil_depth(k)*(1+air_in.w)*air_in.m_dot/(surface_area*viscosity);%Re_a should be between 400 and 1500 % what is delta_coil on page 800? coil thickness?
            D_hdr = 4*coil.min_airflow_area(k)*coil.coil_depth(k)/surface_area;
            C1 = 0.159*(coil.fin_thickness(k)/D_hdr)^(-.065)*(coil.fin_thickness(k)/fin_height)^.141;
            C2 = -0.323*(coil.fin_spacing(k)/fin_height)^.049*(coil.fin_diameter(k)/coil.tube_spacing(k))^0.549*(coil.fin_thickness(k)/coil.fin_spacing(k))^(-0.028);
            f_o = C1*Re_a^C2*air_in.m_dot*Cp_air*Pr^(-2/3);% from referenced paper (doesn't normalize by min_airflow_area like in user manual
%             f_o = C1*Re_a^C2*(air_in.m_dot/coil.min_airflow_area(k))*Cp_air*Pr^(-2/3);
            f_o_w = f_o*(1.425 - 5.1e-4*Re_a + 2.63e-7*Re_a^2);
            n_o = surface_efficiency(f_o,coil.tube_od(k),surface_area,coil.fin_diameter(k),coil.fin_area(k),coil.fin_conductivity(k),coil.fin_thickness(k));
            n_o_w = surface_efficiency(f_o_w,coil.tube_od(k),surface_area,coil.fin_diameter(k),coil.fin_area(k),coil.fin_conductivity(k),coil.fin_thickness(k));
            R_o = 1/(f_o*n_o*surface_area);
            R_o_w = Cp_air/bb/(f_o_w*n_o_w*surface_area);
            R_mf = tube_thickness/(coil.tube_conductivity(k)*coil.area_inside(k));% + 5e-2/coil.area_inside(k);%fouling factor from page 800 is 5e-2 m^2*K/W
            %% iterate mass flow of water to get to desired outlet temperature
%             m_w = m_w_rated;
            Tw_out = max(coil.water_outlet_temperature(k),Tw_in + 1);
            m_w = min(m_w_rated,max(0.1*m_w,sensible_load/(Cp_water*(Tw_out - Tw_in))));
            error = 1;
            count = 0;
            while abs(error)>tolerance
                [air_out.T,air_out.w] = detailed_coil_heat_transfer(m_w,air_in,Tw_in,pipe_area,coil.area_outside(k),coil.tube_id(k),aa,bb,R_mf,R_o,R_o_w);
                dT_dm = (air_out.T - detailed_coil_heat_transfer(m_w*(1+1e-2),air_in,Tw_in,pipe_area,coil.area_outside(k),coil.tube_id(k),aa,bb,R_mf,R_o,R_o_w))/(1e-2*m_w);
                error = air_out.T - T_supply;
                count = count + 1;
                if (error>0 && m_w_rated==m_w) || count>10%exits the loop when it reaches max water flow
                    break
                end
                m_w = min(m_w_rated,max(0.1*m_w,real(m_w + .5*error/dT_dm)));
            end
            air_out.h = psychometric(air_out,'h');
        case 'Water'
            if coil.rated_water_flow(k)>0
                Tw_in0 = coil.water_inlet_temperature(k);
                Ta_in0 = coil.air_inlet_temperature(k);
                Ta_in = air_in.T;
                air_density = 1.204; %kg/m^3
                m_a_rated = coil.rated_air_flow(k)*air_density;
                m_a = air_in.m_dot;
                den_water = 995;%kg/m^3
                m_w_rated = coil.rated_water_flow(k)*den_water;%convert m^3/s to kg/s
                %determine if it has a design capacity from autosizing or just a nominal total capacity
                if isfield(coil,'design_capacity') && length(coil.design_capacity)>=k && coil.design_capacity(k)>0
                    capacity = coil.design_capacity(k);
                else
                    capacity = coil.capacity(k);
                end
                %part of UA calculations depending only on inputs
                if isfield(coil,'air_water_convection_ratio') && ~isnan(coil.air_water_convection_ratio(k))
                    r = coil.air_water_convection_ratio(k);
                    hA_w0 = coil.UA(k)*((r+1)/r);
                    hA_a0 = capacity/(coil.air_outlet_temperature(k) - coil.air_inlet_temperature(k));
                    n_f = r*hA_w0/hA_a0;
                else% need to find root where r and n_f satisfy set of equations
                    switch mode
                        case 'heating'
                            hA_a0 = capacity/(coil.air_outlet_temperature(k) - coil.air_inlet_temperature(k));
                            hA_w0 = capacity/(coil.water_inlet_temperature(k) - coil.water_outlet_temperature(k));
                        case 'cooling'
                            hA_a0 = capacity/(coil.air_inlet_temperature(k) - coil.air_outlet_temperature(k));
                            hA_w0 = capacity/(coil.water_outlet_temperature(k) - coil.water_inlet_temperature(k));
                    end
                    n_f = 1/(hA_a0/coil.UA(k)- hA_a0/hA_w0);%fin efficiency  
                    r = n_f*hA_a0/hA_w0;%air water convection ratio
                end
                x_a = 1 + 4.769e-3*(Ta_in - Ta_in0);
                x_w = 1 + (.014/(1 + .014*Tw_in0))*(Tw_in - Tw_in0);
                if isempty(plant_nodes) 
                    m_w = 0.5*m_w_rated;%start with minimum air flow and half water flow
                    error = 1;
                    while abs(error)>tolerance
                        air_out.T = coil_heat_transfer(x_w,m_w,m_w_rated,hA_w0,x_a,m_a,m_a_rated,hA_a0,n_f,Cp_air,Ta_in,Tw_in);
                        air_out.h = psychometric(air_out,'h');
                        actual_load = air_out.m_dot*(air_out.h - air_in.h);%heat added in W
                        error = (sensible_load - actual_load)/max([sensible_load,actual_load,1]);
                        if error>0 && m_w_rated==m_w%exits the loop when it reaches max water flow
                            break
                        end
                        m_w = min(m_w_rated,m_w*(1+error));
                    end
                    %% If increasing water flow was insufficient, solve for air flow to meet load
                    while abs(error)>tolerance
                        air_out.T = coil_heat_transfer(x_w,m_w,m_w_rated,hA_w0,x_a,air_out.m_dot,m_a_rated,hA_a0,n_f,Cp_air,Ta_in,Tw_in);
                        air_out.h = psychometric(air_out,'h');
                        actual_load = air_out.m_dot*(air_out.h - air_in.h);%heat added in W
                        error = (sensible_load - actual_load)/max([sensible_load,actual_load,1]);
                        if error>0 && m_a_rated==air_out.m_dot%exits the loop when it reaches max water flow
                            break
                        end
                        air_out.m_dot = min(m_a_rated,air_out.m_dot*(1+error));
                    end
                else%if m_w == 0 %first guess solving for m_w
                    %%loop to solve for water flow that gives correct sensible load (page 960 of the manual)
                    %%water flow rate and water outlet temperature change
                    switch mode
                        case 'heating'
                            Tw_out = min(coil.water_outlet_temperature(k),Tw_in - 1);
                            m_w = min(m_w_rated,max(0.1*m_w,sensible_load/(Cp_water*(Tw_in - Tw_out))));
                        case 'cooling'
                            Tw_out = max(coil.water_outlet_temperature(k),Tw_in + 1);
                            m_w = min(m_w_rated,max(0.1*m_w,sensible_load/(Cp_water*(Tw_out - Tw_in))));
                    end
                    error = 1;
                    while abs(error)>tolerance
                        air_out.T = coil_heat_transfer(x_w,m_w,m_w_rated,hA_w0,x_a,m_a,m_a_rated,hA_a0,n_f,Cp_air,Ta_in,Tw_in);
                        switch mode
                            case 'heating'
                                dT_dm = (coil_heat_transfer(x_w,m_w*(1+1e-2),m_w_rated,hA_w0,x_a,m_a,m_a_rated,hA_a0,n_f,Cp_air,Ta_in,Tw_in) - air_out.T)/(1e-2*m_w);
                                error = T_supply - air_out.T;
                            case 'cooling'
                                dT_dm = (air_out.T - coil_heat_transfer(x_w,m_w*(1+1e-2),m_w_rated,hA_w0,x_a,m_a,m_a_rated,hA_a0,n_f,Cp_air,Ta_in,Tw_in))/(1e-2*m_w);
                                error = air_out.T - T_supply;
                        end
                        if error>0 && m_w_rated==m_w%exits the loop when it reaches max water flow
                            break
                        end
                        m_w = min(m_w_rated,max(0.1*m_w,real(m_w + error/dT_dm)));
                    end
                    w_sat = psychometric('Twb',air_out.T,'w');
                    if w_sat<air_out.w
                        air_out.w = w_sat;
                    end
                    air_out.h = psychometric(air_out,'h');
                end
            end
    end
else
    m_w = 0;
end
actual_load = -air_out.m_dot.*(air_out.h - air_in.h);%heat added in W
end %Ends function water_coil

function Ta_out = coil_heat_transfer(x_w,m_w,m_w_rated,hA_w0,x_a,m_a,m_a_rated,hA_a0,n_f,Cp_air,Ta_in,Tw_in)
%find new UA factor
Cp_water = 4186; %J/kg*K
hA_a = x_a*((m_a/m_a_rated)^.8)*hA_a0;
hA_w = x_w*((m_w/m_w_rated)^.85)*hA_w0;
UA = 1/(1/hA_w + 1/(n_f*hA_a));

%capacitance
C_air = Cp_air*m_a;
C_w = Cp_water*m_w;
C_min = min(C_air,C_w);
Z = C_min/max(C_air,C_w);

%Find outlet temperatures
NTU = UA/C_min;
effect = 1 - exp((exp(-Z*NTU.^.78)-1)./(Z*NTU.^(-.22)));
Ta_out = Ta_in + effect*C_min*(Tw_in - Ta_in)/C_air;
end%Ends function coil_heat_transfer

function n_o = surface_efficiency(f_o,tube_od,surface_area,fin_diameter,fin_area,fin_conductivity,fin_thickness)
rho = tube_od/fin_diameter;
fai = (fin_diameter - tube_od)/2*(2*f_o/(fin_conductivity*fin_thickness))^.5;
u_e = fai/(1-rho);
u_b = u_e*rho;
n_fin = -2*rho/(fai*(1+rho))*(besseli(1,u_b)*besselk(1,u_e) - besselk(1,u_b)*besseli(1,u_e))/(besseli(0,u_b)*besselk(1,u_e) - besselk(0,u_b)*besseli(1,u_e));
n_fin = 0.75;%override because getting weird answers for n_fin from modified bessel functions
n_o = 1 - (1 - n_fin)*fin_area/surface_area;
end%Ends function surface_efficiency

function [Ta_out,w_out] = detailed_coil_heat_transfer(m_w,air_in,Tw_in,pipe_area,area_inside,tube_id,aa,bb,R_mf,R_o,R_o_w)
den_water = 995;%kg/m^3
Cp_water = 4186; %J/kg*K
Cp_air = 1006 + 1860*air_in.w; %J/kg*K
V_w = m_w/den_water/pipe_area;%velocity in m/s
f_i = 1.429*(1+0.0146*(Tw_in+273))*V_w^(.8)*tube_id.^(-.2);
R_i = 1/(f_i*area_inside);
UA_dry = (1/(R_i + R_mf + R_o));
%% compute exit temperature with all dry
X = 1/(air_in.m_dot*Cp_air);
Y = -1/(m_w*Cp_water);
Z_d = exp(UA_dry*(X + Y));
K1 = (Z_d-1)/(Z_d +Y/X);
% K1 = X*(Z-1)/(X*Z+Y);
% K2 = (X+Y)*Z/(X*Z+Y);
Tw_2 = Tw_in;
Ta_2 = air_in.T - K1*(air_in.T - Tw_2);
%% Check if exit dewpoint temperature is less than surface temperature
Tdp = psychometric(air_in,'Tdp');
T_s_2 = Tw_2 + ((1/UA_dry)-R_o)/(1/UA_dry)*(Ta_2 - Tw_2);
if Tdp<T_s_2
    %% completely dry coil
    Ta_out = Ta_2;
    w_out = air_in.w;
else
    %% not completely dry, iterate to see how much is wet
    x = (air_in.T - Tdp)/(air_in.T - T_s_2); %Portion of coil that is dry
    error = 1;
    X2 = 1/air_in.m_dot;
    Y2 = -bb/(m_w*Cp_water);
    while abs(error)>1e-3
        [T_s_2,Ta_2,Tw_2] = wet_coil_calc(x,aa,bb,R_i,R_mf,R_o,R_o_w,UA_dry,X,Y,X2,Y2,air_in,Tw_in);
        dTs2_dx = (T_s_2 - wet_coil_calc(x+5e-2,aa,bb,R_i,R_mf,R_o,R_o_w,UA_dry,X,Y,X2,Y2,air_in,Tw_in))/5e-2;
        error = T_s_2 - Tdp;
        if (x <1e-4 && error<0) || dTs2_dx<0
            error = 0;
            x = 0;
        end
        x = max(.25*x,min([4*x,1,x + error/dTs2_dx]));
    end
    Q_w = m_w*Cp_water*(Tw_2 - Tw_in);
    Ha_2 = psychometric('Tdb',Ta_2,'w',air_in.w,'h');%condition at end of dry section (equal to dew point)
    Ha3 = Ha_2 - Q_w/air_in.m_dot;
    Ta_out = psychometric('h',Ha3,'Tdb');
    w_out = psychometric('h',Ha3,'Tdb',Ta_out,'w');
end
end%Ends function detailed_coil_heat_transfer

function [T_s_2,Ta_2,Tw_2] = wet_coil_calc(x,aa,bb,R_i,R_mf,R_o,R_o_w,UA_dry,X,Y,X2,Y2,air_in,Tw_in)
Cp_air = 1006 + 1860*air_in.w; %J/kg*K
if x == 0
    UA_wet = (1/bb/(R_i + R_mf + R_o_w));
    Z_w = exp(UA_wet*(X2 + Y2));
%     K3 = (X2+Y2)/(X2*Z_w + Y2);
    K5 = Z_w*(Y2 + X2)/(X2*Z_w + Y2);
    K6 = Y2*(1 - Z_w)/(bb*(X2*Z_w + Y2));
    Tw_2 = K5*Tw_in + K6*(air_in.h - aa);
    Ta_2 = air_in.T;
    R_x = (R_mf + R_i)/R_o_w;
    T_s_2 = 1/(R_x+1)*(Tw_2 + R_x/bb*(air_in.h - aa));
else
    UA_wet = (1-x)*(1/bb/(R_i + R_mf + R_o_w));
    UA_dry_x = x*UA_dry;
    Z_d = exp(UA_dry_x*(X + Y));
    Z_w = exp(UA_wet*(X2 + Y2));
    K1 = (Z_d-1)/(Z_d +Y/X);
    K5 = Z_w*(Y2 + X2)/(X2*Z_w + Y2);
    K6 = Y2*(1 - Z_w)/(bb*(X2*Z_w + Y2));
    Tw_2 = 1/(1 - K1*K6*Cp_air)*(K5*Tw_in + K6*(air_in.h - Cp_air*air_in.T*K1 -aa));
    Ta_2 = air_in.T - K1*(air_in.T - Tw_2);
    T_s_2 = Tw_2 + ((1/UA_dry_x)-R_o)/(1/UA_dry_x)*(Ta_2 - Tw_2);
end
end%Ends function wet_coil_calc