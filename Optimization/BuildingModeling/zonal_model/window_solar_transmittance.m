function [window_gain,S] = window_solar_transmittance(windows,surfaces,vect_2_sun,weather,internal_irrad)
cos_phi_sun = max(0,vect_2_sun*windows.normal')';
n_win = length(windows.name);
%% transmit/reflect for glazing and simple glazing 
transmit_angle_mod = (windows.transmittance(:,1).*cos_phi_sun.^4 + windows.transmittance(:,2).*cos_phi_sun.^3 + windows.transmittance(:,3).*cos_phi_sun.^2 + windows.transmittance(:,4).*cos_phi_sun + windows.transmittance(:,5));
reflect_angle_mod = (windows.reflectance(:,1).*cos_phi_sun.^4 + windows.reflectance(:,2).*cos_phi_sun.^3 + windows.reflectance(:,3).*cos_phi_sun.^2 + windows.reflectance(:,4).*cos_phi_sun + windows.reflectance(:,5));
reflect = windows.normal_reflectance.*(reflect_angle_mod*ones(1,4));
transmit = [windows.normal_transmittance(:,1),windows.normal_transmittance(:,1),windows.normal_transmittance(:,2),windows.normal_transmittance(:,2)].*(transmit_angle_mod*ones(1,4));

visible_absorptance = max(0,1 - transmit - reflect);%%solve equation on page 290 for absorptance, for single layer it is what is not transmitted or reflected
diffuse_absorptance = visible_absorptance;%cant find any mention describing how to get this paramater for simple window model

G = f_index('Glazing',windows.type);
T_incident = transmit(:,1);%net transmittance factor for direct sunlight
T_incident(G) = transmit(G,1).*transmit(G,3)./(1-reflect(G,3).*reflect(G,2));

window_gain = zeros(length(surfaces.name),1);
S = zeros(4,n_win);
for i = 1:1:n_win   
    switch windows.type{i}
        case 'SimpleGlazingSystem'
            S(1,i) = .5*(weather.DNIWm2*cos_phi_sun(i)*visible_absorptance(i,1) + weather.DHIWm2*diffuse_absorptance(i,1) + internal_irrad.visible(i)*(1 - windows.normal_transmittance(i,1) - windows.normal_reflectance(i,1)));
            S(2,i) = S(1,i) + windows.emittance(i,2).*internal_irrad.long_wave(i);%add interior radiation to inside of 1st pane for simpleglazing
        otherwise
            S(1:2,i) = .5*(weather.DNIWm2*cos_phi_sun(i)*visible_absorptance(i,1) + weather.DHIWm2*diffuse_absorptance(i,1) + internal_irrad.visible(i)*windows.normal_transmittance(i,2)*(1 - windows.normal_transmittance(i,1) - windows.normal_reflectance(i,1)));
            S(3,i) = .5*(weather.DNIWm2*cos_phi_sun(i)*transmit(i,1)*visible_absorptance(i,3) + weather.DHIWm2*windows.normal_transmittance(i,1)*diffuse_absorptance(i,3) + internal_irrad.visible(i)*(1 - windows.normal_transmittance(i,2) + windows.normal_reflectance(i,4)));
            S(4,i) = S(3,i) + windows.emittance(i,2).*internal_irrad.long_wave(i);%add interior radiation to inside of 1st pane for simpleglazing
    end
    %solar distribution (page 212 of reference) 
    z = windows.zone_index(i);
    surf = f_index(z,surfaces.zone_index(1:length(surfaces.height))); %no solar gain onto interior furnishings
    f = f_index('Floor',surfaces.type(surf));
    f_absorb = surfaces.absorptance.visible(surf(f),1);
    f_absorb = f_absorb.*surfaces.area(surf(f))/sum(surfaces.area(surf(f))); %for things like corridors with multiple 'floor' surfaces
    normal2zone = weather.DNIWm2*cos_phi_sun(i)*T_incident(i)*windows.area(i);
    window_gain(surf(f)) = window_gain(surf(f)) + normal2zone*f_absorb;
    window_gain(surf) = window_gain(surf) + (weather.DHIWm2*windows.normal_transmittance(i,1)*windows.area(i) + normal2zone*(1-sum(f_absorb)))*surfaces.area(surf)/sum(surfaces.area(surf));%distribue diffuse sunlight to surfaces in zone
end

% %% glazing
% reflect_angle_mod = zeros(length(windows.name),1);
% transmit_angle_mod = zeros(length(windows.name),1);
% transmit_angle_mod(G) = (windows.transmittance(G,1).*cos_phi_sun(G).^4 + windows.transmittance(G,2).*cos_phi_sun(G).^3 + windows.transmittance(G,3).*cos_phi_sun(G).^2 + windows.transmittance(G,4).*cos_phi_sun(G) + windows.transmittance(G,5));
% reflect_angle_mod(G) = (windows.reflectance(G,1).*cos_phi_sun(G).^4 + windows.reflectance(G,2).*cos_phi_sun(G).^3 + windows.reflectance(G,3).*cos_phi_sun(G).^2 + windows.reflectance(G,4).*cos_phi_sun(G) + windows.reflectance(G,5));
% %% simple glazing(page 295)
% transmit_angle_mod(SG) = cos_phi_sun(SG).*(1+(0.768+0.817*windows.solar_heat_gain(SG).^4).*(1-cos_phi_sun(SG).^2).^(3/2));
% f1 = (((2.403*cos_phi_sun(SG) - 6.192).*cos_phi_sun(SG) + 5.625).*cos_phi_sun(SG) - 2.095).*cos_phi_sun(SG) + 1;
% f2 = (((-1.188*cos_phi_sun(SG) + 2.022).*cos_phi_sun(SG) + 0.137).*cos_phi_sun(SG) - 1.720).*cos_phi_sun(SG);
% reflect_angle_mod(SG) = (f1 + f2.*windows.solar_heat_gain(SG).^.5)./(0.7413 - (0.7396*windows.solar_heat_gain(SG).^.5));