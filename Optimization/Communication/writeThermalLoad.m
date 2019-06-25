function writeThermalLoad(varargin)
%this function is used at E-hub to simulate a thermal load with the fans
global FanPortWrite DateSim
RealData = get_data(DateSim,{'Demand','H'});
fwrite(FanPortWrite,num2str(RealData),'char');