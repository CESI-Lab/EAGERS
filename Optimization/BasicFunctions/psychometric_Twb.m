function Twb = psychometric_Twb(varargin)
%(http://www.ce.utexas.edu/prof/Novoselac/classes/ARE383/Handouts/F01_06SI.pdf)
if length(varargin) == 1
    air = varargin{1};
    T = air.T;
    w = air.w;
else
    T = varargin{1};
    w = varargin{2};
end
Tdb_K = T+273.15; %Tdb (Kelvin)
P = 101.325; % atmospheric pressure (kPa)
satP = exp((-5.8002206e3)./Tdb_K + 1.3914993 - 4.8640239e-2*Tdb_K + 4.1764768e-5*Tdb_K.^2 - 1.4452093e-8*Tdb_K.^3 + 6.5459673*log(Tdb_K))/1000; %saturated water vapor pressure ASHRAE 2013 fundamentals eq. 6 in kPa valid for 0 to 200C
W_s = 0.62198*satP/(P-satP);
Twb = (w*(2501+1.805*T)+1.006*T - 2501*W_s)/(4.186*w-2.381*W_s + 1.006);
end%Ends function psychometric