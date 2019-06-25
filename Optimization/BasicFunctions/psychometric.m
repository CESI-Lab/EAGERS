function out = psychometric(varargin)
%% http://www.ce.utexas.edu/prof/Novoselac/classes/ARE383/Handouts/F01_06SI.pdf
%%%Nomenclature
% Tdb = dry bulb temperature in Celcius
% Twb = wet bulb temperature in Celcius
% Tdp = dew point temperature in Celcius
% w = moisture content (kg H20 / kg dry air)
% h enthalpy J/kg*K
% RH relative humidity
% P pressure in Pa
if length(varargin) == 2
    %structed variable
    air = varargin{1};
elseif length(varargin) == 3 %saturated 
    %input property1, value1, property2, value2, desired property
    air.(varargin{1}) = varargin{2};
else
    %input property1, value1, property2, value2, desired property
    air.(varargin{1}) = varargin{2};
    air.(varargin{3}) = varargin{4};
end
if isfield(air,varargin{end})
    air = rmfield(air,varargin{end});
end
if isfield(air,'T')
    if strcmp(varargin{end},'Tdb')
        air = rmfield(air,'T');
    else
        air.Tdb = air.T;
    end
end
if ~isfield(air,'P')
    air.P = 101325; % atmospheric pressure (Pa)
    % P = 101325*(1-2.25577e-5*z)^5.2559;% atmospheric pressure function of height z(Pa)
    % Tdb = 15-0.0065*z;% atmospheric  temperature function of height z(C)
end
P = air.P; % atmospheric pressure (kPa)
if isfield(air,'h') && isfield(air,'w')
    air.Tdb = (air.h/1000 - air.w*2501)./(1.006 + air.w*1.86);%Temperature in C and converting enthalpy to kJ/kg and humidity in massH2O/massAir
elseif isfield(air,'h') && isfield(air,'Tdb')
    air.w = (air.h/1000 - 1.006*air.Tdb)./(2501 + 1.86*air.Tdb);%kJ/kg 
elseif isfield(air,'RH') && isfield(air,'Tdb')
    satP = sat_press(air.Tdb);
    P_H2O = air.RH/100.*satP; % kPa
    air.w = .621945*(P_H2O/(P-P_H2O));
elseif isfield(air,'w') && isfield(air,'Twb')
    satP = sat_press(air.Twb);%very small error introduced using the wet_bulb temperature for saturation pressure
    W_s = 0.62198*satP./(P-satP);
    air.Tdb = (air.Twb.*(4.186*air.w-2.381*W_s + 1.006) +2501*W_s - air.w*2501)./(1.805*air.w+1.006);
elseif isfield(air,'Twb')
    %% assumes saturated
    air.Tdb = air.Twb;
    satP = sat_press(air.Tdb);
    air.w = 0.62198*satP./(P-satP);
elseif isfield(air,'Tdp')
    %% assumes saturated
    air.Tdb = air.Tdp;
    satP = sat_press(air.Tdp);
    air.w = 0.62198*satP./(P-satP);
elseif isfield(air,'h')
    %% assumes saturated
    air.Tdb = 10;
    error = 1;
    while abs(error)>1e-2
        satP = sat_press(air.Tdb);
        air.w = 0.62198*satP./(P-satP);
%         air.Tdb = (air.h/1000 - air.w*2501)./(1.006 + air.w*1.86);
        h_guess = 1000*(1.006*air.Tdb + air.w.*(2501 + 1.86*air.Tdb));
        error = air.h - h_guess;
        air.Tdb = air.Tdb + 0.3*error/(1006+1860*air.w);
    end
end
if ~isfield(air,varargin{end})
    switch varargin{end}
        case 'h'
            air.h = 1000*(1.006*air.Tdb + air.w.*(2501 + 1.86*air.Tdb));
        case 'Twb'
            %not perfect but close (0.26C) % https://journals.ametsoc.org/doi/full/10.1175/JAMC-D-11-0143.1
            satP = sat_press(air.Tdb);
            P_H2O = P*air.w./(0.62198 + air.w);
            air.RH = P_H2O./satP*100;%Relative humidity in percent
            air.Twb = air.Tdb*atan(0.151977*(air.RH + 8.313659).^.5) ...
                + atan(air.Tdb + air.RH) - atan(air.RH - 1.676331)...
                + 0.00391838*(air.RH).^(3/2).*atan(0.023101*air.RH) - 4.686035;
        case 'Tdp'
            air.Tdp = dew_point(air.Tdb,P,air.w);
        case 'RH'
            satP = sat_press(air.Tdb);
            P_H2O = P*air.w./(0.62198 + air.w);
            air.RH = P_H2O./satP*100;%Relative humidity in percent
        case 'T'
            air.T = air.Tdb;
    end
end
out = air.(varargin{end});
end%Ends function psychometric

function Tdp = dew_point(T,P,w)
P_H2O = P.*w./(0.62198 + w)/1000;% kPa 0.62198 is ratio of molar masses
if T<0
    Tdp = 6.09 + 12.608*log(P_H2O) + 0.4959*log(P_H2O).^2;
else
    Tdp = 6.54 + 14.526*log(P_H2O) + 0.7389*log(P_H2O).^2 + 0.09486*log(P_H2O).^3 + 0.4569*(P_H2O).^0.1984; %Dew point from partial pressure of water using ASHRAE 2013 Fundamentals eqn 39 valid from 0C to 93C
end
end%Ends function dew_point

function satP = sat_press(T)
Tdb_K = T+273.15; %Tdb (Kelvin)
if T<0
    satP = exp((-5.6745359e3)./Tdb_K + 6.3925247 - 9.677843e-3*Tdb_K + 6.22157e-7*Tdb_K.^2 + 2.0747825e-9*Tdb_K.^3 -9.484024e-13*Tdb_K.^4 + 4.1635019*log(Tdb_K)); % Pa, saturated water vapor pressure ASHRAE 2013 fundamentals eq. 6 in kPa valid for 0 to 200C
else
    satP = exp((-5.8002206e3)./Tdb_K + 1.3914993 - 4.8640239e-2*Tdb_K + 4.1764768e-5*Tdb_K.^2 - 1.4452093e-8*Tdb_K.^3 + 6.5459673*log(Tdb_K)); % Pa saturated water vapor pressure ASHRAE 2013 fundamentals eq. 6 in kPa valid for 0 to 200C
end
end%Ends function sat_press