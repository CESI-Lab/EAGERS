function [sunrise,sunset,azimuth,zenith,hour_angle,declination] = solar_calc(long,lat,time_zone,date)
%SOLAR_CALC Calculate position of sun and sunrise and sunset times.
%   [sunrise,sunset,azimuth,zenith] = solar_calc(long,lat,time_zone,date)
%   Calculated using NOAA solar calculations available at:
%   https://www.esrl.noaa.gov/gmd/grad/solcalc/NOAA_Solar_Calculations_day.xls
%   LONG is longitude where + is to the east.
%   LAT is latitude where + is to the north.
%   TIMEZONE is the time zone where + is to the east.
%   DATE is the date in datenum format, i.e. Jan 1 2017 = 736696.
%   SUNRISE and SUNSET are given in fraction of the day, i.e. 6am = 6/24.
%   AZIMUTH and ZENITH are given in degrees.

tpm = date - floor(date);%time past midnight
jd = date + 1721058.5 - time_zone/24;%julian day
jc = (jd-2451545)/36525;
geom_mean_long_sun = mod(280.46646+jc.*(36000.76983+jc*.0003032),360); % degrees
geom_mean_anom = 357.52911+jc.*(35999.05029 - 0.0001537*jc); % degrees
eccent_earth_orbit = 0.016708634-jc.*(0.000042037+0.0000001267*jc);
sun_eq_of_centr = sin(geom_mean_anom/57.2958) .* (1.914602-jc.*(0.004817+0.000014*jc)) + sin(2*geom_mean_anom/57.2958) .* (0.019993-0.000101*jc) + sin(3*geom_mean_anom/57.2958) * 0.000289;
sun_true_long = geom_mean_long_sun + sun_eq_of_centr; % degrees
mean_oblique_corr = 23 + (26 + ((21.448 - jc .* (46.815 + jc .*(.00059-jc*0.001813)))) / 60) / 60; % degrees
sun_app_long = sun_true_long - 0.00569 - 0.00478 * sin((125.04-1934.136*jc)/57.2958); % degrees
oblique_corr = mean_oblique_corr + 0.00256*cos((125.04-1934.136*jc)/57.2958);% degrees
declination = asin(sin(oblique_corr/57.2958) .* sin(sun_app_long/57.2958)); % Solar declination in radians
var_y = tan(oblique_corr/2/57.2958).*tan(oblique_corr/2/57.2958);
eq_of_time = 4 * 360 / (2*pi) * ...
    (var_y .* sin(2*geom_mean_long_sun/57.2958) - 2 * eccent_earth_orbit .* ...
    sin(geom_mean_anom/57.2958) + 4 * eccent_earth_orbit .* var_y .* ...
    sin(geom_mean_anom/57.2958) .* cos(2*geom_mean_long_sun/57.2958) - .5 * ...
    var_y .* var_y .* sin(4*geom_mean_long_sun/57.2958) - 1.25 * ...
    eccent_earth_orbit .* eccent_earth_orbit .* sin(2*geom_mean_anom/57.2958));
ha_sunrise = 360 / (2*pi) * acos(cos(90.833/57.2958) ./ ...
    (cos(lat/57.2958).*cos(declination)) - tan(lat/57.2958) .* ...
    tan(declination)); % degrees (sunlight hours)
solar_noon = (720-4*long - eq_of_time + time_zone*60) / 1440;
sunrise = solar_noon - ha_sunrise * 4 / 1440; % Local Sidereal Time (LST)
sunset = solar_noon + ha_sunrise * 4 / 1440; % Local Sidereal Time (LST)

tst = mod(tpm*1440 + eq_of_time + 4*long - 60*time_zone, 1440); % True Solar Time in min
hour_angle = (tst/4+180)/57.2958; % radians
hour_angle(tst>0) = hour_angle-2*pi;% radians
zenith = 360 / (2*pi) * acos(sin(lat/57.2958) .* sin(declination) +  cos(lat/57.2958) .* cos(declination) .* cos(hour_angle)); % degrees
% ElevationAngle = 90 - Zenith;
% Refraction = zeros(length(ElevationAngle),1);
% Refracion(ElevationAngle<85 & ElevationAngle>5) = (58.1./tan(ElevationAngle(ElevationAngle<85 & ElevationAngle>5)/57.2958)-0.07./tan(ElevationAngle(ElevationAngle<85 & ElevationAngle>5)/57.2958).^3 +.000086./tan(ElevationAngle(ElevationAngle<85 & ElevationAngle>5)/57.2958).^5)/3600;
% Refracion(ElevationAngle<=5 & ElevationAngle>-.575) = (1735+(ElevationAngle(ElevationAngle<=5 & ElevationAngle>-.575).*(-518.2+ElevationAngle(ElevationAngle<=5 & ElevationAngle>-.575).*(103.4+ElevationAngle(ElevationAngle<=5 & ElevationAngle>-.575).*(-12.79+ElevationAngle(ElevationAngle<=5 & ElevationAngle>-.575)*0.711)))))/3600;
% Refracion(ElevationAngle<=.575) = -20.722./tan(ElevationAngle(ElevationAngle<=.575)/57.2958)/3600;
% CorrectedElevationAngle = ElevationAngle+Refraction;

ang = 360/(2*pi)*real(acos(((sin(lat/57.2958).*cos(zenith/57.2958)) - sin(declination))./(cos(lat/57.2958).*sin(zenith/57.2958))));
azimuth = mod(540 - ang,360); % degrees clockwise from N
azimuth(hour_angle>0) = mod(ang(hour_angle>0)+180,360); % degrees clockwise from N
end % Ends function solar_calc
