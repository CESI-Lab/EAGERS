function [T_z,T_s,m_v,T_w] = update_model_states_implicit(building,T,gains,window_gain,dt,w,mix,plenum,infiltration,central_flow,direct_flow,weather)
%% Analytic equation soln (eqn 2.11 in engineering reference) for zone temperature
%% Finite diference Implicit solution for surface temperatures
air_density = 1.204; %kg/m^3 %standard air density (1.204 kg/m3) adjusted for the local barometric pressure (standard barometric pressure corrected for altitude, ASHRAE 1997 HOF pg. 6.1).
cp_zone = air_density*(1006 + 1860*w.zone); %J/m^3*K
cp_window = cp_zone(building.windows.zone_index);
z = length(T.zone);
c_t = 1*building.zones.volume;%sensible heat multiplier
s_area = building.surfaces.area;
%%m*Cp & h*A*T terms for grouping all temperature and humidity states as (T.zone;m_v_zone] for analytic solution  
return_flow = sum(mix+plenum,2) + infiltration + central_flow + direct_flow;
den_w = air_density*return_flow;
den_w(den_w==0) = air_density*building.zones.volume(den_w==0)/dt;%adjustment for zones with zero air flow of any kind into them
cap = [cp_zone.*c_t;air_density*building.zones.volume];%capacitance [zone_thermal;zone_moisture];
m_cp_no_mix = cp_zone.*(infiltration.*T.air_zone +  central_flow.*T.central + direct_flow.*T.direct);%flows in m^3/s converted to kg/s*J/(kg*K)*K = W 
h2o_no_mix = gains.zone_latent'/2264705 + air_density*(infiltration.*w.air +  central_flow.*w.central + direct_flow.*w.direct);%flows in m^3/s converted to kg/s, everything into kg/s of water 

%determine how small of time steps to take
sub_steps = 1;
tc = dt/sub_steps;
errorTz = 0;
errorTs = 0;
errorW = 0;
Ts_j0 = T.surf;
for i = 1:1:sub_steps %loop through the wall at a faster time freq, so that internal temperatures distribute
    %initial guess for next step: interpolate from T.zone and T.zone_est from previous call to update_model_states, and add error in sub_steps
    g = i/sub_steps;
    Tz_i = (1-g)*T.zone + g*T.zone_est;%zone temperature at next time step
    Ts_i = (1-g)*T.surf + g*T.surf_est;
    mv_i = (1-g)*w.zone + g*w.zone_est;
    Tz_j1 = Tz_i + .5*errorTz;
    Ts_j1 = Ts_i + .5*errorTs;
    mv_j1 = mv_i + .5*errorW;
    %% converge to a solution for the sub-step
    max_e = 100*building.site.temp_tol;
    count = 0;
    old_Ts = Ts_j1;
    while max_e>0.1*building.site.temp_tol
        %surface temperature portion
        [Ts_j1b,h_interior,max_e,count,old_Ts] = surface_implicit(building,weather,T,gains,window_gain,Tz_j1,Ts_j1,Ts_j0,tc,count,old_Ts);
        %zone temperature and humidity portion
        h_windows = window_convection(building.windows,(T.windows_int+273),(Tz_j1(building.windows.zone_index)+273),cp_window);
        h_A = building.ctf.map_sur2zone*(h_interior.*s_area) +  building.ctf.map_win2zone*(h_windows.*building.windows.area);
        h_A_Twin = building.ctf.map_win2zone*(h_windows.*building.windows.area.*T.windows_int);
        den_z = h_A + cp_zone.*return_flow;% convection + inter-zone mixing + infiltration + HVAC
        den = [den_z;den_w];
        m_cp_T = cp_zone.*((mix+plenum)*Tz_j1) + m_cp_no_mix;
        h_A_T = building.ctf.map_sur2zone*(h_interior.*s_area.*Ts_j1b(building.ctf.sur_state1)) + h_A_Twin;
        term1_z = gains.zone_sensible' + h_A_T + m_cp_T;% internal gains + convection*T.surface + inter-zone mixing*T.source + infiltration*T.air + HVAC*T.supply_zone   
        term1_w = air_density*((mix+plenum)*mv_j1) + h2o_no_mix;%flows in m^3/s converted to kg/s, everything into kg/s of water 
        term1 = [term1_z;term1_w];
        new = ([Tz_j1;mv_j1] - term1./den).*exp(-den./cap*tc) + term1./den;%Analytic equation soln (eqn 2.11 in engineering reference)
        max_e = max([max(abs(Ts_j1-Ts_j1b)),max(abs(new(1:z)-Tz_j1)),100*max(abs(mv_j1 - new(z+1:end)))]);
        Tz_j1 = new(1:z);%zone temperature at end of sub-step
        Ts_j1 = .25*Ts_j1+.75*Ts_j1b;
        mv_j1 = new(z+1:end);%zone humidity at end of sub-step
    end
%     %% check energy conservation
%     sig = 5.67e-8;% W/m^2*K^4 Stephan-Boltzman const
%     n_sur = length(building.ctf.sur_state1);
%     n_win = length(T.windows_int);
%     T_surf_interior = Ts_j1(building.ctf.sur_state1);
%     T_surf_exterior = Ts_j1(building.ctf.sur_state2);
%     dT = T_surf_interior - building.ctf.map_sur2zone'*Tz_j1;%interior surface convection = h*A*(Tsur - Tzone) 
%     q_radiative = sig*building.zones.view_factors(1:n_sur,1:n_sur).*(building.ctf.interior_absorb*ones(1,n_sur)).*((T_surf_interior+273).^4 - (T_surf_interior+273)'.^4);%heat transfer from i to j: h_ij*(Ti - Tj) = Q_i -->j
%     q_rad2windows = sig*building.zones.view_factors(1:n_sur,n_sur+1:end).*(building.ctf.interior_absorb*ones(1,n_win)).*((T_surf_interior+273).^4 - (T.windows_int+273)'.^4);%heat transfer from i to j: h_ij*(Ti - Tj) = Q_i -->j     
%     q_interior = (gains.interior_surface' + window_gain)./building.surfaces.area - h_interior.*dT - sum(q_radiative,2) - sum(q_rad2windows,2);
%     [h_exterior,h_sky] = exterior_convection(building,Ts_j1(building.ctf.sur_state2),T.air_surface,weather,T.sky,building.convection.exterior);
%     q_exterior = gains.exterior_surface' + h_exterior.*(T.exterior - T_surf_exterior) + h_sky.*(T.sky - T_surf_exterior) + building.ctf.ground_cond.*(T.ground - T_surf_exterior);
%     B = building.ctf.capacitance./building.ctf.surf_area/tc;
%     error_energy = (-building.ctf.A*Ts_j1 + building.ctf.interior_surface*q_interior + building.ctf.exterior_surface*q_exterior) - (Ts_j1 - Ts_j0).*B;% (Energy into node - dT*capacitance ) = error in energy balance (W)
%     error_temp_change_frac = error_energy./B./(Ts_j1 - Ts_j0);
%     error_energy = error_energy.*building.ctf.surf_area;
%     %%%---%%
    
    Ts_j0 = Ts_j1;%used as 'previous surface temperature for next sub-step
    errorTz = Tz_j1-Tz_i;
    errorTs = Ts_j1-Ts_i;
    errorW = mv_j1-mv_i;
%%need fix so that when raining outside exterior surface temperature = T wet bulb (page 91 of manual)
end
T_z = new(1:z);
T_s = Ts_j1;
m_v = new(z+1:end);
T_w = T.windows_int;
end%Ends function update_model_states

function [Ts_j1,h_interior,max_e,count,old_Ts] = surface_implicit(building,weather,T,gains,window_gain,Tz_j1,Ts_j1,Ts_j0,tc,count,old_Ts)
%surface temperature portion
sig = 5.67e-8;% W/m^2*K^4 Stephan-Boltzman const
max_e = 100*building.site.temp_tol;
n_sur = length(building.ctf.sur_state1);
n_win = length(T.windows_int);
%scale implicit matrix by time constant
n = length(building.ctf.surf_area);
B = building.ctf.capacitance./building.ctf.surf_area/tc;
A = building.ctf.A + spdiags(B,0,n,n);

r = .25;
if count == 0
    change = 8;
else
    change = min(4,max(1e-2,abs(Ts_j1 - Ts_j0)));
end
while max_e>3*building.site.temp_tol
        T_surf_interior = Ts_j1(building.ctf.sur_state1);
        T_surf_exterior = Ts_j1(building.ctf.sur_state2);
        dT = T_surf_interior - building.ctf.map_sur2zone'*Tz_j1;%interior surface convection = h*A*(Tsur - Tzone) 
        h_interior = natural_convection(building.surfaces.normal,dT,[],building.convection.interior);
        q_radiative = sig*building.zones.view_factors(1:n_sur,1:n_sur).*(building.ctf.interior_absorb*ones(1,n_sur)).*((T_surf_interior+273).^4 - (T_surf_interior+273)'.^4);%heat transfer from i to j: h_ij*(Ti - Tj) = Q_i -->j
        q_rad2windows = sig*building.zones.view_factors(1:n_sur,n_sur+1:end).*(building.ctf.interior_absorb*ones(1,n_win)).*((T_surf_interior+273).^4 - (T.windows_int+273)'.^4);%heat transfer from i to j: h_ij*(Ti - Tj) = Q_i -->j     
        q_interior = (gains.interior_surface' + window_gain)./building.surfaces.area - h_interior.*dT - sum(q_radiative,2) - sum(q_rad2windows,2);
        [h_exterior,h_sky] = exterior_convection(building,Ts_j1(building.ctf.sur_state2),T.air_surface,weather,T.sky,building.convection.exterior);
        q_exterior = gains.exterior_surface' + h_exterior.*(T.exterior - T_surf_exterior) + h_sky.*(T.sky - T_surf_exterior) + building.ctf.ground_cond.*(T.ground - T_surf_exterior);
        b = B.*Ts_j0 + building.ctf.interior_surface*q_interior + building.ctf.exterior_surface*q_exterior;
        new_Ts = A\b;
        error = new_Ts - Ts_j1;
        old_error = old_Ts- Ts_j1;
        max_e  = max(abs(error));
        ds = r*error;
        overshoot = sign(error)~=sign(old_error);
        ds(overshoot) = .33*ds(overshoot);
        ds = max(-.5*change,min(.5*change,ds));
        old_Ts = new_Ts;
        Ts_j1 = Ts_j1 + ds;
        change = min(4,max(1e-2,abs(Ts_j1 - Ts_j0)));
        count = count + 1;
end
% %%checking things
% error_energy = (-building.ctf.A*Ts_j1 + building.ctf.interior_surface*q_interior + building.ctf.exterior_surface*q_exterior) - (Ts_j1 - Ts_j0).*B;% (Energy into node - dT*capacitance ) = error in energy balance (W)
% error_temp_change_frac = error_energy./B./(Ts_j1 - Ts_j0);
% error_energy = error_energy.*building.ctf.surf_area;
end%Ends function surface_implicit