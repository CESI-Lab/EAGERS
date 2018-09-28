function h = psychometric_h(varargin)
if length(varargin) == 1
    air = varargin{1};
    T = air.T;
    w = air.w;
else
    T = varargin{1};
    w = varargin{2};
end
h = (1.006*T + w.*(2501 + 1.86*T));%kJ/kg %ASHRAE 2013 fundamentals eq. 32 (http://www.ce.utexas.edu/prof/Novoselac/classes/ARE383/Handouts/F01_06SI.pdf)
end%Ends function psychometric