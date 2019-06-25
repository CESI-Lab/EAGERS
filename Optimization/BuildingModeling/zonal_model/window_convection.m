function h_i = window_convection(windows,T_int,T_i,cp_air)
%% correlations for changing tilt angle (page 362)
%% note temperatures must be in Kelvin!!
air_density = 1.204; %kg/m^3 %standard air density (1.204 kg/m3) adjusted for the local barometric pressure (standard barometric pressure corrected for altitude, ASHRAE 1997 HOF pg. 6.1).
g = 9.81; %m/s^2 gravity
deltaT = T_int-T_i;
nusselt = zeros(length(windows.name),1);
tilt = abs(atand((windows.normal(:,1).^2+windows.normal(:,2).^2).^.5./windows.normal(:,3)));
tilt(deltaT>0) = 180 - tilt(deltaT>0);
T_mf = T_i + 0.25*deltaT;
mu = 3.723e-6 + (4.94e-8)*T_mf;
lambda = 2.873e-3 + (7.76e-5)*T_mf;
Ra_h = air_density^2*windows.height.^3*g.*cp_air.*abs(deltaT)./(T_mf.*mu.*lambda);
Ra_cv = 2.5e5*(exp(.72*tilt)./sind(tilt)).^(1/5);

r1 = tilt<15;
nusselt(r1) = 0.13*Ra_h(r1).^(1/3);
r2 = tilt>=15 & tilt<=90 & Ra_h<=Ra_cv;
nusselt(r2) = 0.56*(Ra_h(r2).*sind(tilt(r2))).^(1/4);
r3 = tilt>=15 & tilt<=90 & Ra_h>Ra_cv;
nusselt(r3) = 0.13*(Ra_h(r3).^(1/3) - Ra_cv(r3).^(1/3)) + 0.56*(Ra_cv(r3).*sind(tilt(r3))).^(1/4);
r4 = tilt>90 & tilt<=179;
nusselt(r4) = 0.56*(max(min(Ra_h(r4),1e11),1e5).*sind(tilt(r4))).^(1/4);
r5 = tilt>179;
nusselt(r5) = 0.58*(min(Ra_h(r5),1e11)).^(1/5);

h_i = nusselt.*(lambda./windows.height);
end%Ends function window_convection