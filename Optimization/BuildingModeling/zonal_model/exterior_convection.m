function [h_s,h_sky] = exterior_convection(building,T_exterior_sur,Tair,weather,T_sky,method)
sig = 5.67e-8;% W/m^2*K^4 Stephan-Boltzman const
exterior = building.surfaces.exterior;
F_grnd = 0.5*(1-exterior.cos_phi);
F_sky = 0.5*(1+exterior.cos_phi);
a = f_index('Ground',exterior.boundary);
F_sky(a) = 0;
F_grnd(a) = 0;
beta = F_sky.^.5;
h_n = natural_convection(exterior.normal,T_exterior_sur - Tair,[],method);
leward = leward_calc(weather.Wdirdegrees,exterior.normal);
wind_speed = wind_speed_calc(weather.Wspdms,building.ctf.s_height,building.site.terrain{1});
switch method
    case 'DOE-2'
        h_glass = (h_n.^2 + (3.26*wind_speed.^.89).^2).^.5;
        h_glass(leward) = (h_n(leward).^2 + (3.55*wind_speed(leward).^.617).^2).^.5;
        h_wind = h_n + exterior.roughness_factor.*(h_glass-h_n);
    case {'Detailed';'BLAST';'TARP'}
        h_wind = h_n + 2.537*(1-.5*leward).*exterior.roughness_factor.*(exterior.perimeter.*wind_speed./exterior.area).^.5;
end
h_wind(a) = 0;
h_ground = exterior.thermal_absorptance*sig.*F_grnd.*((T_exterior_sur+273).^4 - (Tair+273).^4)./(T_exterior_sur - Tair);
h_air = exterior.thermal_absorptance*sig.*F_sky.*(1-beta).*((T_exterior_sur+273).^4 - (Tair+273).^4)./(T_exterior_sur - Tair);
h_sky = exterior.thermal_absorptance*sig.*F_sky.*beta.*((T_exterior_sur+273).^4 - (T_sky+273).^4)./(T_exterior_sur - T_sky);
alt = (T_exterior_sur == Tair);
h_ground(alt) = exterior.thermal_absorptance(alt)*sig.*F_grnd(alt).*4.*(Tair(alt)).^3;%avoid nan by using derivative
h_air(alt) = exterior.thermal_absorptance(alt)*sig.*F_sky(alt).*(1-beta(alt)).*4.*(Tair(alt)).^3;%avoid nan by using derivative
alt = (T_exterior_sur == T_sky);
h_sky(alt) = exterior.thermal_absorptance(alt)*sig.*F_sky(alt).*beta(alt).*4.*(T_sky).^3;%avoid nan by using derivative

h_s = h_wind + h_air + h_ground;
end%Ends function exterior convection