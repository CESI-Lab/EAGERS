function solution = solver_nrel(gen,options,ic,T0,forecast,scale_cost)
dt = options.Resolution*3600;% time step for discrete model, seconds
n_s = options.Horizon/options.Resolution;% # of time steps in optimization
N_sch = options.Horizon; % assumes hourly scheduling
forecast.internal_gain = forecast.Demand.IntGain;

% Initial conditions for building and equipment status
Pfc_f0 = ic(1)*1000;%initial fuel cell power converted from kW to W
Pfc_h0 = 0;
Eb0 = ic(5)*1000*3600;%initial state of charge converted from kWh to J
u0 = 0;

dt_h = 60*60; % Time step for scheduling period
N_perhour = dt_h/dt; % Number of thermal dynamics model steps per hour
N_state = 3600/dt*N_sch; % Number of thermal dynamics model steps


Toa = forecast.Weather.Tdb(N_perhour:N_perhour:n_s); % outside air temperature
Pl_h = forecast.internal_gain(N_perhour:N_perhour:n_s)*1e3;%internal heating gains in watts
Pl_e = forecast.Demand.E(N_perhour:N_perhour:n_s)*1e3;%non HVAC electric load in watts

% alpha = 0.4; % outdoor air ratio
m_oa = 0.3; % outdoor air flow rate, kg/s


% % Building zone parameters
zone_param.Cp = 1e3; % specific heat of air
zone_param.Cr = 5e6;
zone_param.Cw = 3.75e8;
zone_param.R1 = 0.031e-3;
zone_param.R2 = 0.5898e-3;
zone_param.Ts = 12.8; % supply air temperature
zone_param.T_approx = 22; % approximation for room temperature for complexity

COP = gen(4).Output.Cooling(end); % chiller COP
eta_fc_e = gen(3).Output.Electricity(end); % electrical efficiency of the FC
Pfc_f_max = gen(3).Size*1000/eta_fc_e; % FC fuel consumption capacity


DischCurrent = gen(5).VariableStruct.PeakDisch.*gen(5).Size/gen(5).VariableStruct.Voltage*1000;
DischResistScaled = (100/DischCurrent)*gen(5).VariableStruct.DischResist*(1/1000); %Scale so the loss of power is equivelant to that specified at 100Amps
DischVoltLoss = DischCurrent*DischResistScaled; %keep in mind when calculating loss as function of discharge current
ChargeCurrent = gen(5).VariableStruct.PeakCharge*gen(5).Size/gen(5).VariableStruct.Voltage*1000;
ChargeResistScaled = (100/ChargeCurrent)*gen(5).VariableStruct.ChargeResist*(1/1000); %Scale so the loss of power is equivelant to that specified at 100Amps
ChargeVoltLoss = ChargeCurrent*ChargeResistScaled;
eta_c = gen(5).VariableStruct.Voltage/(gen(5).VariableStruct.Voltage+ChargeVoltLoss); %charging efficiency
eta_d = (gen(5).VariableStruct.Voltage-DischVoltLoss)/gen(5).VariableStruct.Voltage; %discharging efficiency
Eb_max = gen(5).Size*(gen(5).VariableStruct.MaxDOD/100)*1000*3600; % usable battery capacity (J)
Pb_max = (DischCurrent*gen(5).VariableStruct.Voltage*eta_d); % battery ramping capacity (J/s)
eta_fc_h = 0.35; % heat efficiecny of the FC
eta_boil = 0.9; % boiler efficiency
a_f = 2.225e+03; % fan power coefficient
b_f = -2.977e+03;

kas = 1e-4; % coefficent to calculate temperature buffer from ancillary service
c_up = 1; % FC start-up cost
c_down = 1; % FC shut-down cost

T_high = 23; % temperature bounds
T_low = 21;
Tref_high = 32; % temperature set-points bounds
Tref_low = 12;
% J = zeros(nS,1);

Pfc_f_min = 0.2*Pfc_f_max; % FC fuel consumption minimum operating point

m_min = 2.23; % supply air flow rate limits
m_max = 10;
FC_ramp_up = 0.25*Pfc_f_max; % FC ramp rate limits
FC_ramp_down = 0.25*Pfc_f_max;
min_up = 6; % FC minimum up time
min_down = 6; % FC minimum down time

[num_vec,den_vec]=zone_model(zone_param,dt); % Model parameter generator

%% EDC optimization
cvx_begin

    % Decision variables:
    % zone temperatue set-point, battery charge, battery discharge, supply
    % air flow rate, FC fuel consumption set-point, FC heat output
    % set-point, total reheat, AS capacity
    variables T_ref(N_sch) Pb_in(N_sch) Pb_out(N_sch) m(N_sch) Pfc_f_ref(N_sch) Pfc_h_ref(N_sch) Prh(N_sch) Ras(N_sch)
    % State variables:
    % zone temperature, battery SOC
    variables T(N_state) Eb(N_state)
    % Interger variables:
    % FC on/off, FC start up, FC shut down
    variable u(N_sch+max(min_up,min_down)) binary
    variable u_up(N_sch+max(min_up,min_down)) binary
    variable u_down(N_sch+max(min_up,min_down)) binary

% calculate cost
for i = 1:N_sch
    T_ref_i = T_ref(i);
    Pb_in_i = Pb_in(i);
    Pb_out_i = Pb_out(i);
    m_i = m(i);
    Pfc_f_ref_i = Pfc_f_ref(i);
    Pfc_h_ref_i = Pfc_h_ref(i);
    Prh_i = Prh(i);
    Ras_i = Ras(i);
    u_i = u(i);
    u_up_i = u_up(i);
    u_down_i = u_down(i);
    
% Total cost incurred from various sources (electricity, natural gas, fuel cell, AS)
    % Electricity cost
    Pcc_i = (m_i-m_oa)*zone_param.Cp*(zone_param.T_approx-zone_param.Ts) + m_oa*zone_param.Cp*(Toa(i)-zone_param.Ts); % cooling coil heat exchange
    Pch_i = 1/COP*Pcc_i; % chiller power

    Pf_i = a_f*m_i + b_f; % fan power
    
    Pfc_e_i = Pfc_f_ref_i*eta_fc_e; % FC elec output
    
    Pl_e_i = Pl_e(i); % non-HVAC elec load
    
    P_e_i = Pch_i + Pf_i + Pb_in_i + Pl_e_i- Pfc_e_i - Pb_out_i; % grid power
    
    J_e_i = scale_cost(i*N_perhour,1)*P_e_i;
    % Gas cost
    P_ng_i = Pfc_f_ref_i + (Prh_i-Pfc_h_ref_i)/eta_boil;
    J_gas_i = scale_cost(i*N_perhour,2)*P_ng_i;    
    % AS payment
    J_as_i = -scale_cost(i*N_perhour,3)*Ras_i;    
    % FC on/off
    J_fc_com_i = c_up*u_up_i + c_down*u_down_i;   
    % Total cost
    J(i,1) = J_e_i + J_gas_i + J_as_i + J_fc_com_i;
    
end


minimize ( sum(abs(J)) )

subject to

for i_sch = 1:N_sch  % (N_sch=24)
     T_ref_i = T_ref(i_sch);
    Pb_in_i = Pb_in(i_sch);
    Pb_out_i = Pb_out(i_sch);
    m_i = m(i_sch);
    Pfc_f_ref_i = Pfc_f_ref(i_sch);
    Pfc_h_ref_i = Pfc_h_ref(i_sch);
    Prh_i = Prh(i_sch);
    Ras_i = Ras(i_sch);
    u_i = u(i_sch);

    if i_sch == 1
        Pfc_f_ref_iminus = Pfc_f0;
        Pfc_h_ref_iminus = Pfc_h0;
        u_iminus = u0;
    else
        Pfc_f_ref_iminus = Pfc_f_ref(i_sch-1);
        Pfc_h_ref_iminus = Pfc_h_ref(i_sch-1);
        u_iminus = u(i_sch-1);
    end

    
    % Limits
    0 <= Pfc_h_ref_i <= Pfc_f_ref_i*eta_fc_h;
    0 <= Prh_i - Pfc_h_ref_i;
    0 <= Ras_i <= 1e3;
    m_min <= m_i <= m_max;
    Pfc_f_min*u_i <= Pfc_f_ref_i <= Pfc_f_max*u_i;
    Pfc_f_min*1 <= Pfc_f_ref_i <= Pfc_f_max*1;
    -FC_ramp_down <= Pfc_f_ref_i - Pfc_f_ref_iminus <= FC_ramp_up;
    0 <= Pb_in_i <= Pb_max;
    0 <= Pb_out_i <= Pb_max;
    Tref_low <= T_ref_i <= Tref_high;
    
    % FC commitment
    for i_up = 1:min_up % FC minimum up time
        -u_iminus + u_i - u(i_sch+i_up) <= 0;
    end
    for i_down = 1:min_down % FC minimum down time
        u_iminus - u_i + u(i_sch+i_down) <= 1;
    end

    -u_iminus + u_i - u_up(i_sch) <= 0; % FC start up
    u_iminus - u_i - u_down(i_sch) <= 0; % FC shut down
            
    % States in faster time step
    for i_state = 1:N_perhour
        T_i = T((i_sch-1)*N_perhour+i_state);
        Eb_i = Eb((i_sch-1)*N_perhour+i_state);

        if i_sch == 1 && i_state == 1
            T_iminus = T0;
            T_iminus2 = T0;
            Eb_iminus = Eb0;
        elseif i_sch == 1 && i_state == 2
            T_iminus = T(1);
            T_iminus2 = T0;
            Eb_iminus = Eb((i_sch-1)*N_perhour+i_state-1);
        else
            T_iminus = T((i_sch-1)*N_perhour+i_state-1);
            T_iminus2 = T((i_sch-1)*N_perhour+i_state-2);
            Eb_iminus = Eb((i_sch-1)*N_perhour+i_state-1);    
        end
        Pl_h_iminus = Pl_h(i_sch);
        Pl_h_iminus2 = Pl_h(i_sch);
        Toa_iminus = Toa(i_sch);
        Toa_iminus2 = Toa(i_sch);

        T_new = -den_vec(3,2)*T_iminus - den_vec(3,3)*T_iminus2 + num_vec(1,2)*m_i + num_vec(2,2)*(Pl_h_iminus+Prh_i) + num_vec(3,2)*Toa_iminus + num_vec(1,3)*m_i + num_vec(2,3)*(Pl_h_iminus2+Prh_i) + num_vec(3,3)*Toa_iminus2;
        Eb_new = Eb_iminus + dt*eta_c*Pb_in_i - dt*1/eta_d*Pb_out_i;
        
        T_new == T_i;
                         
        Eb_new == Eb_i
        
        T_low + kas*Ras_i <= T_i <= T_high - kas*Ras_i;
        0 <= Eb_i <= Eb_max;
        
    end
    
    T_new == T_ref_i; % temperature reach set-point
    
    
end

Eb_new == Eb0; % final battery SOC

cvx_end

solution.Dispatch = zeros(n_s+1,5); 
solution.heat_recovery = zeros(n_s,1); 
solution.Buildings.Heating = zeros(n_s,1); 
solution.Dispatch(1,:) = ic;
solution.Dispatch(2:end,5)= Eb/1000/3600;   %battery state of charge (converted from J to kWh)
for t = 1:1:N_sch
    index = 2+(t-1)*N_perhour:1+t*N_perhour;
    solution.Dispatch(index,3)= Pfc_f_ref(t)*eta_fc_e/1000; % FC elec output; % Fuel cell electrical output
    solution.Dispatch(index,4)= ((m(t)-m_oa)*zone_param.Cp*(zone_param.T_approx-zone_param.Ts) + m_oa*zone_param.Cp*(Toa(t)-zone_param.Ts))/1000;   % Chiller output
    solution.Dispatch(index,1) = (1000/COP*solution.Dispatch(index,4) + (a_f*m(t) + b_f) + Pb_in(t) - Pb_out(t) + Pl_e(t) - Pfc_f_ref(t)*eta_fc_e + (Prh(t)-Pfc_h_ref(t))/eta_boil - Ras(t))/1000;%power to/from the electric grid in kW  
    solution.Buildings.Heating(index,1) = Prh(t);
    solution.heat_recovery(index,1)  = Pfc_h_ref(t);
end
solution.Buildings.Cooling = solution.Dispatch(:,4);
solution.Buildings.Temperature = T;
end%Ends function solver_nrel

function [num_vec,den_vec]= zone_model(zone_param,dt)
% Formulate MISO model for zone temperature
% inputs: u = [m, (Plh+Prh), Toa]
A = zeros(2,2);
A(1,1) = 1-dt/(zone_param.Cr*zone_param.R1);
A(1,2) = dt/(zone_param.Cr*zone_param.R1);
A(2,1) = dt/(zone_param.Cw*zone_param.R1);
A(2,2) = 1-dt/(zone_param.Cw*zone_param.R1)-dt/(zone_param.Cw*zone_param.R2);

B = zeros(2,3);
B(1,1) = dt/zone_param.Cr*zone_param.Cp*(zone_param.Ts-zone_param.T_approx);
B(1,2) = dt/zone_param.Cr;
B(2,3) = dt/(zone_param.Cw*zone_param.R2);

C = [1,0];

D = zeros(1,3);

num_vec = zeros(3,3);
den_vec = zeros(3,3);
for i_ss = 1:3
    [num,den] = ss2tf(A,B,C,D,i_ss);
    num_vec(i_ss,:)=num;
    den_vec(i_ss,:)=den;
end
end%ends function zone_modle