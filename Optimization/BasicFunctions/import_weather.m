function weather = import_weather(filename)
weather.filename = filename;
%% Initialize variables.
delimiter = ',';
%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
[filepath,name,ext] = fileparts(filename);
weather.Timestamp = linspace(datenum([2017,1,1,1,0,0]),datenum([2018,1,1]),8760)';
if strcmpi(ext,'.csv')
    formatSpec = '%s%s%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%s%f%f%s%f%f%s%f%f%s%f%f%s%f%f%s%f%f%s%f%f%s%f%f%s%f%f%s%f%f%s%f%f%s%f%f%s%f%f%f%s%f%f%s%f%[^\n\r]';
    startRow = 3;
    endRow = inf;
    dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines', startRow(1)-1, 'ReturnOnError', false);
    
    Headers = {'DateMMDDYYYY';'TimeHHMM';'ETRWm2';'ETRNWm2';'GHIWm2';'GHIsource';'GHIuncert';'DNIWm2';'DNIsource';'DNIuncert';'DHIWm2';'DHIsource';'DHIuncert';'GHillumlx';'GHillumsource';'Globalillumuncert';'DNillumlx';...
        'DNillumsource';'DNillumuncert';'DHillumlx';'DHillumsource';'DHillumuncert';'Zenithlumcdm2';'Zenithlumsource';'Zenithlumuncert';'TotCldtenths';'TotCldsource';'TotClduncertcode';'OpqCldtenths';'OpqCldsource';...
        'OpqClduncertcode';'DrybulbC';'Drybulbsource';'Drybulbuncertcode';'DewpointC';'Dewpointsource';'Dewpointuncertcode';'RHum';'RHumsource';'RHumuncertcode';'Pressurembar';'Pressuresource';'Pressureuncertcode';...
        'Wdirdegrees';'Wdirsource';'Wdiruncertcode';'Wspdms';'Wspdsource';'Wspduncertcode';'Hvism';'Hvissource';'Hvisuncertcode';'CeilHgtm';'CeilHgtsource';'CeilHgtuncertcode';'Pwatcm';'Pwatsource';'Pwatuncertcode';...
        'AODunitless';'AODsource';'AODuncertcode';'Albunitless';'Albsource';'Albuncertcode';'Lprecipdepthmm';'Lprecipquantityhr';'Lprecipsource';'Lprecipuncertcode';'PresWthMETARcode';'PresWthsource';'PresWthuncertcode';};
    % H_keep = [5,8,11,14,17,20,23,26,29,32,35,38,41,44,47,50,53,56,59,62,65,66];
    H_keep = [5,8,11,32,35,38,41,44,47];
    for i = 1:1:length(H_keep)
        weather.(Headers{H_keep(i)}) = dataArray{H_keep(i)};
    end
    weather.PressurePa = 100*weather.Pressurembar;
elseif strcmpi(ext,'.epw')
    formatSpec = '%f%f%f%f%f%s%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';%
    startRow = 9;
    endRow = inf;
    dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines', startRow(1)-1, 'ReturnOnError', false);
    dataArray(6:end-1) = dataArray(7:end);%eliminate column A1 with uncertainties
    Headers = {'Year';'Month';'Day';'Hour';'Minute';'DrybulbC';'DewpointC';'RHum';'PressurePa';'ETRWm2';'ETRNWm2';'GHinfaredWm2';'GHIWm2';'DNIWm2';'DHIWm2';'GHillumlx';'DNillumlx';'DHillumlx';...
        'Zenithlumcdm2';'Wdirdegrees';'Wspdms';'TotCldtenths';'OpqCldtenths';'Hvism';'CeilHgtm';};
    H_keep = [6,7,8,9,13,14,15,20,21,22,23];
    for i = 1:1:length(H_keep)
        weather.(Headers{H_keep(i)}) = dataArray{H_keep(i)};
    end
    weather.TotCldtenths = weather.TotCldtenths/10;
    weather.OpqCldtenths = weather.OpqCldtenths/10;
end
%% Close the text file.
fclose(fileID);
%% Field descriptions for TMY3 data
% 1 Date MM/DD/YYYY -- Date of data record 
% 2 Time HH:MM -- Time of data record (local standard time) 
% 3 Hourly extraterrestrial radiation on a horizontal surface Watt-hour per square meter 1 Wh/m2 Amount of solar radiation received on a horizontal surface at the top of the atmosphere during the 60-minute period ending at the timestamp 
% 4 Hourly extraterrestrial radiation normal to the sun Watt-hour per square meter 1 Wh/m2 Amount of solar radiation received on a surface normal to the sun at the top of the atmosphere during the 60- minute period ending at the timestamp 
% 5 Global horizontal irradiance Watt-hour per square meter 1 Wh/m2 Total amount of direct and diffuse solar radiation received on a horizontal surface during the 60-minute period ending at the timestamp 
% 6 Global horizontal irradiance source flag 1-2 -- See Table 1-4 
% 7 Global horizontal irradiance uncertainty Percent 1% Uncertainty based on random and bias error estimates – see NSRDB User’s Manual (Wilcox, 2007b) 
% 8 Direct normal irradiance Watt-hour per square meter 1 Wh/m2 Amount of solar radiation (modeled) received in a collimated beam on a surface normal to the sun during the 60- minute period ending at the timestamp 
% 9 Direct normal irradiance source flag 1-2 -- See table 1-4 
% 10 Direct normal irradiance uncertainty Percent 1% Uncertainty based on random and bias error estimates – see NSRDB User’s Manual (Wilcox, 2007b) 
% 11 Diffuse horizontal irradiance Watt-hour per square meter 1 Wh/m2 Amount of solar radiation received from the sky (excluding the solar disk) on a horizontal surface during the 60-minute period ending at the timestamp 
% 12 Diffuse horizontal irradiance source flag 1-2 -- See Table 1-4 
% 13 Diffuse horizontal irradiance uncertainty Percent 1% Uncertainty based on random and bias error estimates – see NSRDB User’s Manual (Wilcox, 2007b) 4 
% 14 Global horizontal illuminance Lux 100 lx Average total amount of direct and diffuse illuminance received on a horizontal surface during the 60-minute period ending at the timestamp 
% 15 Global horizontal illuminance source flag 1-2 -- See Table 1-4 
% 16 Global horizontal illuminance uncertainty Percent 1% Uncertainty based on random and bias error estimates – see section 2.10) 
% 17 Direct normal illuminance Lux 100 lx Average amount of direct normal illuminance received within a 5.7° field of view centered on the sun during 60- minute period ending at the timestamp 
% 18 Direct normal illuminance source flag 1-2 -- See Table 1-4 
% 19 Direct normal illuminance uncertainty Percent 1% Uncertainty based on random and bias error estimates – see section 2.10) 
% 20 Diffuse horizontal illuminance Lux 100 lx Average amount of illuminance received from the sky (excluding the solar disk) on a horizontal surface during the 60-minute period ending at the timestamp 
% 21 Diffuse horizontal illuminance source flag 1-2 -- See Table 1-4
% 22 Diffuse horizontal illuminance uncertainty Percent 1% Uncertainty based on random and bias error estimates – see section 2.10) 
% 23 Zenith luminance Candela per square meter 10 cd/m2 Average amount of luminance at the sky's zenith during the 60- minute period ending at the timestamp 
% 24 Zenith luminance source flag 1-2 -- See Table 1-4 
% 25 Zenith luminance uncertainty Percent 1% Uncertainty based on random and bias error estimates – see section 2.10) 
% 26 Total sky cover Tenths of sky 1 tenth Amount of sky dome covered by clouds or obscuring phenomena at the time indicated 
% 27 Total sky cover flag (source) See Table 1-5 
% 28 Total sky cover flag (uncertainty) See Table 1-6 
% 29 Opaque sky cover Tenths of sky 1 tenth Amount of sky dome covered by clouds or obscuring phenomena that prevent observing the sky or higher cloud layers at the time indicated 5 
% 30 Opaque sky cover flag (source) See Table 1-5 
% 31 Opaque sky cover flag (uncertainty) See Table 1-6 
% 32 Dry-bulb temperature Degree C 0.1° Dry-bulb temperature at the time indicated 
% 33 Dry-bulb temperature flag (source) See Table 1-5 
% 34 Dry-bulb temperature flag (uncertainty) See Table 1-6 
% 35 Dew-point temperature Degree C 0.1° Dew-point temperature at the time indicated 
% 36 Dew-point temperature flag (source) See Table 1-5 
% 37 Dew-point temperature flag (uncertainty) See Table 1-6 
% 38 Relative humidity Percent 1% Relative humidity at the time indicated
% 39 Relative humidity flag (source) See Table 1-5 
% 40 Relative humidity flag (uncertainty) See Table 1-6 
% 41 Station pressure Millibar 1 mbar Station pressure at the time indicated 
% 42 Station pressure flag (source) See Table 1-5 
% 43 Station pressure flag (uncertainty) See Table 1-6 
% 44 Wind direction Degrees from north (360° = north; 0° = undefined, calm) 10° Wind direction at the time indicated 
% 45 Wind direction flag (source) See Table 1-5 
% 46 Wind direction flag (uncertainty) See Table 1-6 
% 47 Wind speed Meter/second 0.1 m/s Wind speed at the time indicated 
% 48 Wind speed flag (source) See Table 1-5 
% 49 Wind speed flag (uncertainty) See Table 1-6 
% 50 Horizontal visibility Meter* 1 m Distance to discernable remote objects at the time indicated (7777 = unlimited) 
% 51 Horizontal visibility flag (source) See Table 1-5 
% 52 Horizontal visibility flag (uncertainty) 6 
% 53 Ceiling height Meter* 1 m Height of the cloud base above local terrain (77777 = unlimited)
% 54 Ceiling height flag (source) See Table 1-5 
% 55 Ceiling height flag (uncertainty) See Table 1-6 
% 56 Precipitable water Centimeter 0.1 cm The total precipitable water contained in a column of unit cross section extending from the earth's surface to the top of the atmosphere 
% 57 Precipitable water flag (source) See Table 1-5 
% 58 Precipitable water flag (uncertainty) See Table 1-6 
% 59 Aerosol optical depth, broadband [unitless] 0.001 The broadband aerosol optical depth per unit of air mass due to extinction by the aerosol component of the atmosphere 
% 60 Aerosol optical depth, broadband flag (source) See Table 1-5 
% 61 Aerosol optical depth, broadband flag (uncertainty) See Table 1-6 
% 62 Albedo [unitless] 0.01 The ratio of reflected solar irradiance to global horizontal irradiance 
% 63 Albedo flag (source) See Table 1-5 
% 64 Albedo flag (uncertainty) See Table 1-6 
% 65 Liquid precipitation depth Millimeter* 1 mm The amount of liquid precipitation observed at the indicated time for the period indicated in the liquid precipitation quantity field 
% 66 Liquid precipitation quantity Hour* 1 hr The period of accumulation for the liquid precipitation depth field
% 67 Liquid precipitation depth flag (source) See Table 1-5 
% 68 Liquid precipitation depth flag (uncertainty) See Table 1-6

%% Fields for epw format

% 1: Year
% 2: Month
% 3: Day
% 4: Hour
% 5: Minute
% A1: Data Source and Uncertainty Flags
% \note Initial day of weather file is checked by EnergyPlus for validity (as shown below)
% \note Each field is checked for "missing" as shown below. Reasonable values, calculated
% \note values or the last "good" value is substituted.
% 6: Dry Bulb Temperature
% units: C
% minimum: > -70
% maximum: < 70
% missing: 99.9
% 7: Dew Point Temperature
% units: C
% minimum: > -70
% maximum: < 70
% missing: 99.9
% 8: Relative Humidity
% missing: 999.
% minimum: 0
% maximum: 110
% 9: Atmospheric Station Pressure
% units: Pa
% missing: 999999.
% minimum: > 31000
% maximum: < 120000
% 10: Extraterrestrial Horizontal Radiation
% units: Wh/m2
% missing: 9999.
% minimum: 0
% 11: Extraterrestrial Direct Normal Radiation
% units: Wh/m2
% missing: 9999.
% minimum: 0
% 12: Horizontal Infrared Radiation Intensity
% units: Wh/m2
% missing: 9999.
% minimum: 0
% 13: Global Horizontal Radiation
% units: Wh/m2
% missing: 9999.
% minimum: 0
% 14: Direct Normal Radiation
% units: Wh/m2
% missing: 9999.
% minimum: 0
% 15: Diffuse Horizontal Radiation
% units: Wh/m2
% missing: 9999.
% minimum: 0
% 16: Global Horizontal Illuminance
% units: lux
% missing: 999999.
% minimum: 0
% 17: Direct Normal Illuminance
% units: lux
% missing: 999999.
% minimum: 0
% 18: Diffuse Horizontal Illuminance
% units: lux
% missing: 999999.
% minimum: 0
% 19: Zenith Luminance
% units: Cd/m2
% missing: 9999.
% minimum: 0
% 20: Wind Direction
% units: degrees
% missing: 999.
% minimum: 0
% maximum: 360
% 21: Wind Speed
% units: m/s
% missing: 999.
% minimum: 0
% maximum: 40
% 22: Total Sky Cover
% missing: 99
% minimum: 0
% maximum: 10
% 23: Opaque Sky Cover (used if Horizontal IR Intensity missing)
% missing: 99
% minimum: 0
% maximum: 10
% 24: Visibility
% units: km
% missing: 9999
% 25: Ceiling Height
% units: m
% missing: 99999
% 26: Present Weather Observation
% 27: Present Weather Codes
% 28: Precipitable Water
% units: mm
% missing: 999
% 29: Aerosol Optical Depth
% units: thousandths
% missing: .999
% 30: Snow Depth
% units: cm
% missing: 999
% 31: Days Since Last Snowfall
% missing: 99
% 32: Albedo
% missing: 999
% 33: Liquid Precipitation Depth
% units: mm
% missing: 999
% 34 Liquid Precipitation Quantity
% units: hr
% missing: 99