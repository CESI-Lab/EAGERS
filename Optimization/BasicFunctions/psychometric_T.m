function T = psychometric_T(varargin)
if length(varargin) == 1
    air = varargin{1};
    h = air.h;
    w = air.w;
else
    h = varargin{1};
    w = varargin{2};
end
T = (h - w*2501)./(1.006 + w*1.86);%Temperature in C using enthalpy in kJ/kg and humidity in massH2O/massAir
end%Ends function psychometric