function solution = solver_nrel(gen,building,options,T0,forecast,scale_cost)

% Load inputs for optimization
%   Electricity price, natural gas price, ancillary service rate, outdoor
%   air temp, Q internal, solar irradiance
c_e = scale_cost(:,1); % electricity price
c_ng = scale_cost(:,1); % natural gas price
if isfield(forecast,'ancillary_service')%placeholder, we dont currently have any datasets to load into here
    c_as = forecast.ancillary_service; % ancillary service rate
else
    c_as = 0; % ancillary service rate
end
T_amb = forecast.Weather.DrybulbC; % outdoor air temperature
Sol_rad = forecast.Weather.DNIWm2;% solar irradiance
Qint = forecast.Building.InternalGains;  % Q internal

% Load building parameters -- HARD CODE THIS PARAMETERS
%   ROM, fan, chiller, boiler, supply air flow limits, reheat limits,
%   ancillary service buffer coefficent, outdoor air ratio
A = building.A; % ROM 
B = building.B;
E = building.E;
Fan_base = building.Fan_base; % 2nd-order poly, 3 coefficients
Fan_bot = building.Fan_bot;
Fan_mid = building.Fan_mid;
Fan_top = building.Fan_top;
Chiller = building.Chiller; % 3 coefficients
Boiler = building.Boiler; % 3 coefficients
m_min = building.m_min; % supply air flow limits
m_max = building.m_max;
Prh_min = building.Prh_min; % reheat limits
Prh_max = building.Prh_max;
kas = building.kas; % ancillary service buffer coefficient 
alpha_oa =building.alpha_oa; % outdoor air ratio

% Load sizable equipment parameters
%   Fuel cell: capacity (min/max power), electric efficiency, heat efficiency, ramp rate limits, start-up/shut-down costs, minimum up/down time
%   Battery: capacity (power and energy), charging/discharging efficiency
%   Thermal storage: capacity (power and energy), charging/discharging
%   efficiency, static loss
Pfc_f_max = gen(3).Size*gen(3).Output.Capacity(end)/gen(3).Output.DirectCurrent(end); % FC fuel consumption capacity
Pfc_f_min = gen(3).Size*gen(3).Output.Capacity(2)/gen(3).Output.DirectCurrent(2); % FC fuel consumption minimum operating point
eta_fc_e = mean(gen(3).Output.DirectCurrent(2)); % electrical efficiency of the FC
eta_fc_h = mean(gen(3).Output.Heat(2)); % heat efficiecny of the FC
FC_ramp_up = gen(3).VariableStruct.dX_dt; % FC ramp rate limits
FC_ramp_down = gen(3).VariableStruct.dX_dt;
c_up = gen(3).VariableStruct.StartCost; % FC start-up cost

%% we don't have the following in our generator models (They won't show up in GUI for editing and other generators in library would not be compatible)
c_down = gen(3).VariableStruct.StopCost; % FC shut-down cost
min_up = gen(3).VariableStruct.MinUp; % FC min up time
min_down = gen(3).VariableStruct.MinDown; % FC mim down time
%%--%%

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

Pth_max = gen(7).VariableStruct.DischRatePerc/100*gen(7).Size; % thermal storage power capacity
Eth_max = gen(7).Size; % thermal storage energy capacity
eta_th_d = gen(7).VariableStruct.DischargeEff; % thermal storage discharging efficiency
eta_th_c = gen(7).VariableStruct.ChargeEff; % thermal storage charging efficiency
eta_th_loss = gen(7).VariableStruct.SelfDischarge; % thermal storage static loss

% Load optimization setup parameters
%   Temperature comfortable range, temp setpoint range, time step, number
%   of scheduling period
T_high = forecast.Building.Tmin; % temperature bounds
T_low = forecast.Building.Tmax;
Tref_high = forecast.Building.Tmin; % temperature set-points bounds
Tref_low = forecast.Building.Tmax;
Ts = forecast.Building.cooling_air_temperature; % supply air temperature
T_approx = T0;% approximation for room temperature for complexity
dt = options.Resolution*3600; % Time step for discrete thermal dynamics model, seconds
N_sch = options.Horizon;%(Assuming hourly) % Number of scheduling period

% Load initial condition
%   Zone temperature, FC power and heat output, batter and thermal storage
%   SOC
X0 = T0;   % Initial temperature values
Pfc_f0 = gen(1).CurrentState(1)*1000; % FC fuel consumption 
Pfc_h0 = gen(1).CurrentState(1)*1000*eta_fc_h/eta_fc_e; % FC heat output
Eb0 = gen(5).CurrentState(1)*1000*3600; % battery SOC
Eth0 = gen(7).CurrentState(1)*1000*3600; % thermal storage SOC


%% EDC optimization

cvx_begin

cvx_solver sedumi

% Decision variables:
% zone temperatue set-point, battery charge, battery discharge, supply
% air flow rate, FC fuel consumption set-point, FC heat output
% set-point, total reheat, AS capacity
variables u(N_sch,5) m_z(N_sch,19) Prh(N_sch,16) T_ref(N_sch,16) 
% u =[Pb_in,Pb_out,Pfc_f_ref,Pfc_h_ref,Ras]
% m_z: Zone airflow rate, '19' indicates 19 zones, similar for reheat (Prh)
% T_ref: setpoint temperature (No plenums (3) included, hence '16')

% State variables:
% zone temperature, battery SOC
variables x(N_sch,19,N_perhour) Eb(N_sch,1,N_perhour)  % Three dimensional vector (second column represents 'size')
%('19' includes zones (Basement, Bottomfloor(6),Middlefloor(6),Topfloor(6))

expression Pcc_total(N_sch)
expression Ppre_total(N_sch)
expression Pf_total(N_sch)
expression Pcc_i(3)
expression J_i(N_sch)

% calculate cost
for i = 1:N_sch
    
    % Total cost incurred from various sources (electricity, natural gas, fuel cell, AS)
    % Electricity cost
    
  %% Chiller cost
    Tmix_i = Tmix_vec(i);
    for n_f = 1:N_AHU-1
        if Tmix_i >= Zoneparam.Ts(n_f)
            Pcc_i(n_f,1) = sum(m_z(i,6*n_f-4:6*n_f))*Zoneparam.Cp*(Tmix_i-Zoneparam.Ts(n_f,1));
        else
            Pcc_i(n_f,1) = 0;
        end
    end
     
   
        if Tmix_i >= Zoneparam.Ts(n_f)
            Pcc_b = m_z(i,1)*Zoneparam.Cp*(Tmix_i-Zoneparam.Ts(n_f,1));
        else
            Pcc_b = 0;
        end
  
    
    Pcc_total(i,1) = sum(Pcc_i) + Pcc_b;   % addition of the CAV system chiller power consumption
    
    Pcc_chiller = Chiller.d0+Chiller.d1*Tamb(i)+Chiller.d2*(Pcc_total(i,1));  % Regression equation for chiller
    
    Pch_i = Pcc_chiller; % chiller power

    %%  Fan
    for n_f = 1:N_AHU
        if n_f ==1
            Pf_i(n_f,1) = C_f(1,n_f)+C_f(2,n_f)*m_z(i,n_f);
        else
        Pf_i(n_f,1) = C_f(1,n_f)+C_f(2,n_f)*sum(m_z(i,6*(n_f-1)-4:6*(n_f-1))) + C_f(3,n_f)*power(sum(m_z(i,6*(n_f-1)-4:6*(n_f-1))),2);%m_s(i,n_f).^2;  % m_i is the supply air flow rate 
        end
    end
    
    Pf_total(i,1) = sum(Pf_i);
    
    %
    Pfc_e_i = u(i,3)*eta_fc_e; % FC elec output
    
    P_e_i = Pch_i + Pf_total(i,1) + u(i,1) - Pfc_e_i - u(i,2); % grid power
    
    J_e_i = c_e(i)*P_e_i;
    
    Ppre_total(i,1) = sum(Prh(i,:),2);
    
    Pcc_boiler(i,1) = Boiler.e0+Boiler.e1*Tamb(i)+Boiler.e2*(Ppre_total(i,1)-u(i,4));
    
    P_ng_i = u(i,3) + Pcc_boiler(i,1);
    J_gas_i = c_ng(i)*P_ng_i;
    % AS payment
    J_as_i = -c_as(i)*u(i,5);
    % Total cost
    J(i,1) = J_e_i + J_gas_i + J_as_i;
    
end


minimize ( sum(J) )


subject to

for i_sch = 1:N_sch  % (nS=24)
    
    
    T_ref_i = T_ref(i_sch,:);
    m_i = m_z(i_sch,:);
    Prh_i = Prh(i_sch,:);
    Pb_in_i = u(i_sch,1);
    Pb_out_i = u(i_sch,2);
    Pfc_f_ref_i = u(i_sch,3);
    Pfc_h_ref_i = u(i_sch,4);
    Ras_i = u(i_sch,5);
    
    if i_sch == 1
        Pfc_f_ref_iminus = Pfc_f0;
        Pfc_h_ref_iminus = Pfc_h0;
    else
        Pfc_f_ref_iminus = u((i_sch-1),3);
        Pfc_h_ref_iminus = u((i_sch-1),4);
    end
    
    %Limits
    0 <= Pfc_h_ref_i <= Pfc_f_ref_i*eta_fc_h;
    0 <= sum(Prh(i_sch,:)) - Pfc_h_ref_i; 
    0 <= Ras_i <= 0;
    m_min <= m_i' <= m_max;
    Pfc_f_min*1 <= Pfc_f_ref_i <= Pfc_f_max*1;
    -FC_ramp_down <= Pfc_f_ref_i - Pfc_f_ref_iminus <= FC_ramp_up;
    0 <= Pb_in_i <= Pb_max;
    0 <= Pb_out_i <= Pb_max;
    Tref_low <= T_ref_i' <= Tref_high;
    Prh >= 0;
  
    
    % States in faster time step
    
    for i_state = 1:N_perhour
        
        T_i = x(i_sch,:,i_state)';
        
        Eb_i = Eb(i_sch,1,i_state);
        
        if i_sch == 1 && i_state == 1
            X = X0;
            Eb_iminus = Eb0;
        elseif i_state == 1
            X = x(i_sch-1,:,N_perhour)';
            Eb_iminus = Eb(i_sch-1,1,N_perhour);
        else
            X = x(i_sch,:,i_state-1)';
            Eb_iminus = Eb(i_sch,1,i_state-1);
        end
        
        U = [T_ref_i'; m_i';Prh_i'];  % size(U)=51=[16;19;16]
        D = [Tamb(i_sch);Sol_rad(i_sch);Qint(i_sch,:)';1]; %size(D)=19=[1;1;16;1]
        
   
        Xnew=A*X+B*U+E*D; % X=[Basement,Bot_floor_zones,Bot_plenum,Mid_floor_zones,Mid_plenum,Top_floor_zones,Top_plenum] total = 19 states
        Eb_new = Eb_iminus + dt*eta_c*Pb_in_i - dt*1/eta_d*Pb_out_i;
        
        Xnew == T_i;
        Eb_new == Eb_i
        
     
        T_low + kas*Ras_i/3 <= [T_i(1:6);T_i(8:12);T_i(14:18)] <= T_high - kas*Ras_i/3;
        0 <= Eb_i <= Eb_max;
        
    end
    
    [Xnew(1:6);Xnew(8:12);Xnew(14:18)]' == T_ref_i; % temperature reach set-point
    
    
end

Eb_new == Eb0; % final battery SOC

cvx_end


%% Post process

% Collect all decision variables as output
%   For all zones: temperature setpoint, supply air flow rate, reheat,
%   ancillary service capacity, actual zone temperature 
%   Operation of FC, battery, thermal storage, SOC of battery and thermal
%   storage

solution.Buildings.Temperature_set = T_ref; % zone temperature setpoint
solution.Buildings.Flow = m_z; % zone supply air flow rate
solution.Buildings.Reheat = Prh; % zone reheat
solution.AncilaryService = u(:,5); % ancillary service capacity
solution.Buildings.Temperature = x; % actual zone temperature
for t = 1:1:N_sch
    index = 2+(t-1)*N_perhour:1+t*N_perhour;
    solution.Dispatch(index,3) = u(t,3)*eta_fc_e/1000; % FC elec output; % Fuel cell electrical output
    solution.heat_recovery(index,1) = u(t,4); % FC heat recovered by building
    solution.Dispatch(index,5) = Eb(t); % battery SOC
    solution.Dispatch(index,7) = Eth(t); % thermal storage SOC
    %% Chiller/Boiler setpoints?
    
    %% Power purcahsed from Grid?
    
end
%% don't need these  (redundant)
u_FC; % FC on/off status %should be infered from fuel consumption
u(:,1); % battery charging
u(:,2); % battery discharging
u(:,6); % thermal storage charging
u(:,7); % thermal storage discharging

