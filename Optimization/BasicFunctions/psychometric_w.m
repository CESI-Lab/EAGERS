function w = psychometric_w(T,h)
%T = dry bulb temperature in C
%h = enthalpy (kJ/kg)
w = (h - 1.006*T)./(2501 + 1.86*T);%kJ/kg %ASHRAE 2013 fundamentals  (http://www.ce.utexas.edu/prof/Novoselac/classes/ARE383/Handouts/F01_06SI.pdf)
end%Ends function psychometric