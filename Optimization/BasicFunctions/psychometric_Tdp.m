function Tdp = psychometric_Tdp(varargin)
%% Ambient dewpoint
%(http://www.ce.utexas.edu/prof/Novoselac/classes/ARE383/Handouts/F01_06SI.pdf)
if length(varargin) == 1
    air = varargin{1};
    T = air.T;
    w = air.w;
else
    T = varargin{1};
    w = varargin{2};
end
P_H2O = 101.325*w*(0.62198 + w);
if T<0
    Tdp = 6.09 + 12.608*log(P_H2O) + 0.4959*log(P_H2O).^2;
else
    Tdp = 6.54 + 14.526*log(P_H2O) + 0.7389*log(P_H2O).^2 + 0.09486*log(P_H2O).^3 + 0.4569*(P_H2O).^0.1984; %Dew point from partial pressure of water using ASHRAE 2013 Fundamentals eqn 39 valid from 0C to 93C
end
end%Ends function psychometric