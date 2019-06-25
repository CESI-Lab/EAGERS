function [T,T_windows_int,Q_windows,window_h] = window_temperature(building,vect_2_sun,T,T_sky,T_zone,T_surf,cp_zone,weather,internal_irrad)
%% all temperatures in Kelvin for this function
% T is a matrix where each row holds the temperature states for a single window. 
% The first column is the inside surface of the innermost pane, second column is outside of first pane, third is inside of 2nd pane...
sig = 5.67e-8;% W/m^2*K^4 Stephan-Boltzman const
windows = building.windows;
n_sur = length(building.surfaces.name);

T_new = T;
SG = f_index('SimpleGlazingSystem',windows.type);
G = f_index('Glazing',windows.type);
k_air = 26.3e-3*windows.gap_thickness;
T_int = T(:,2);
T_int(G) = T(G,4);
T_i = T_zone(windows.zone_index)+273;%Temperature on interior window suface in Kelvin
cp = cp_zone(windows.zone_index);
T_o = weather.DrybulbC - 0.0065*windows.height+273;%outdoor air temp next to each window
T_s = T_surf(building.ctf.sur_state1)+273;%zone interior surface temperatures
[~,S] = window_solar_transmittance(building.windows,building.surfaces,vect_2_sun,weather,internal_irrad);

wind_speed = wind_speed_calc(weather.Wspdms,windows.elevation,building.site.terrain{1});
leward = leward_calc(weather.Wdirdegrees,windows.normal);%orientation relative to wind (surface normal is outward facing)
Z = sig*windows.emittance(:,2).*windows.emittance(:,3)./(1-(1-windows.emittance(:,2)).*(1-windows.emittance(:,3)));
F_sky = 0.5*(1+windows.cos_phi);
Eo = sig.*((0.5*(1-windows.cos_phi)).*T_o.^4 + F_sky.^1.5.*T_sky.^4 + F_sky.*(1-F_sky.^.5).*T_o.^4);
Ei_s = sig.*(building.zones.view_factors(1:n_sur,n_sur+1:end)'*(building.surfaces.area_rad.*T_s.^4))./building.windows.area;%internal radiation from other surfaces to window
vf_w = building.zones.view_factors(n_sur+1:end,n_sur+1:end);%window 2 window view factors
error = 1;
count = 0;
while error>0.1*building.site.temp_tol && count<10
    h_n_ext = natural_convection(windows.normal, T(:,1) - T_o,[],building.convection.exterior);
    switch building.convection.exterior
        case 'DOE-2'
            h_o = (h_n_ext.^2 + ((3.26*wind_speed).^(.89)).^2).^.5;
            h_o(leward) = (h_n_ext(leward).^2 + (3.55*wind_speed(leward).^.617).^2).^.5;
        case {'Detailed';'BLAST';'TARP'}
            h_forced = 2.537*(1-.5*leward).*(windows.perimeter.*wind_speed./windows.area).^.5;
            h_o = h_n_ext + h_forced;
    end
    h_i = window_convection(windows,T_int,T_i,cp);

    for i = 1:1:length(windows.name)
        h_r = windows.emittance(i,:)*sig.*T(i,:).^3;
        k = windows.thermal_conductivity(i,:)./windows.thickness(i,:);
        Ei = Ei_s(i) + sig*(vf_w(i,:)*T_int.^4);
        switch windows.type{i}
            case 'SimpleGlazingSystem'
                e = windows.emittance(i,1:2);
                A = [-h_r(1)-k(1)-h_o(i),       k(1);...
                           k(1)         ,   -h_r(2)-k(1)-h_i(i)];
                B = [-e(1)*Eo(i) - h_o(i)*T_o(i);...
                     -e(2)*Ei - h_i(i)*T_i(i);] - S(1:2,i);
                T_new(i,1:2) = (A\B)';
            otherwise
                if k(2) == 0 || isinf(k(2))
                    disp('single layer window glazing?')
                end
                e = [windows.emittance(i,1),windows.emittance(i,4)];
                h_r(2) = Z(i)*T(i,2)^3;
                h_r(3) = Z(i)*T(i,3)^3;
                h1 = k_air(i,1);
                A = [-h_r(1)-k(1)-h_o(i),       k(1)       ,        0       ,      0;...
                           k(1)         ,   -h_r(2)-k(1)-h1,      h1+h_r(3) ,      0;...
                             0          ,       h1+h_r(2)  , -h1-k(2)-h_r(3),     k(2);...
                             0          ,         0        ,       k(2)     , -h_r(4)-k(2)-h_i(i);];
                B = [-e(1)*Eo(i) - h_o(i)*T_o(i);...
                                 0;...
                                 0;...
                     -e(2)*Ei - h_i(i)*T_i(i);] - S(:,i);
                T_new(i,1:4) = (A\B)';
        end
    end
    error = max(max(abs(T_new - T)));
    r = .4;
    T = (1-r)*T_new + r*T;
    T_int(SG) = T(SG,2);
    T_int(G) = T(G,4);
    count = count+1;
end

window_h = h_i;
Q_windows = building.ctf.map_win2zone*(h_i.*windows.area.*(T_int - T_i));% + output.solar_gain;
T_windows_int = T_int - 273;%convert back to Celcius 
end%Ends function window_temperature