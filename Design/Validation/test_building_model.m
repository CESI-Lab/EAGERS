function [errors,Rsquare] = test_building_model(building,weather,date,mtr)

% building is the EAGERS building structure, weather is a structure of dry
% bulb, wetbuld and relative humididty, and solar irradiation,
% date is the timsteps to be compared (must be sequential).
location = struct('Latitude',40, 'Longitude',-105, 'TimeZone',-7);
dt = date(2:end) - date(1:end-1);
dt = [dt(1);dt(1);dt];
Ti = 17.4;
weather_interp = interpolate_weather(weather,date);
SG = solar_gain(building,date,location,weather_interp);
ExternalGains = SG.Walls + SG.Roof;
    %%warm-up period
    xF = nnz(date<date(1)+7);
    wuDate = date(1:xF);
    wuWeather = interpolate_weather(weather,wuDate);
    wuSG = solar_gain(building,wuDate,location,wuWeather);
    wuLoads = building_loads(building,wuDate,wuSG);
    wuExternalGains = wuSG.Walls + wuSG.Roof;
    [~,~,~,Tzone,Twall,~] = building_profile(building,wuDate,wuLoads.InternalGains,wuExternalGains,wuWeather.Tdb,wuWeather.RH,Ti,Ti);


if ~isempty(mtr)
    Cooling = mtr.CoolingElectricity*building.VariableStruct.COP_C;
    Heating = (mtr.HeatingElectricity + mtr.HeatingGas)*building.VariableStruct.COP_H;
    Damper = ones(length(dt),1);
    Tzone_Eplus = building_simulate(building,weather_interp.Tdb,weather_interp.RH,dt,mtr.ZoneTotalInternalTotalHeatingEnergy,ExternalGains,Cooling,Heating,mtr.AirFlow,Damper,Tzone(end,1),Twall(end,1));
    
    %%plot building model response with EnergyPlus Heating & cooling
    mtr.ZoneMeanAirTemperature
end
%% Run EAGERS method% Compare a forecasted building to it's simulated response
B_loads = building_loads(building,date,SG);    
[Cooling, Heating, Fan_Power,Tzone,Twall,Damper] = building_profile(building,date,B_loads.InternalGains,ExternalGains,weather_interp.Tdb,weather_interp.RH,Tzone(end,1),Twall(end,1));    

%% Simulate response to this Heating/Cooling Profile
AirFlow = Fan_Power/Build.VariableStruct.FanPower;
[Tzone_S,Twall_S] = building_simulate(building,weather_interp.Tdb,weather_interp.RH,dt,B_loads.InternalGains,ExternalGains,Cooling,Heating,AirFlow,Damper,Ti,Ti);
%% Plotting
errors = zeros(9, 1);
Rsquare = zeros(9, 1);
fprintf('Total %% errors:\n')

% %compare equipment
% figTitle = 'Equipment';
% figNum = 1;
% ylab = 'Electric Load (kW)';
% ePlus = mtr.InteriorEquipmentElectricity(xi:xf);
% eagers = B_loads.Equipment;
% [errors(figNum),Rsquare(figNum)] = elec_load_plot_setup(figTitle, figNum, Xe, ePlus, DofY, eagers, dt, dtE,ylab);
% 
% %compare interior lighting
% figTitle = 'Interior Lighting';
% figNum = 2;
% ylab = 'Electric Load (kW)';
% ePlus = mtr.InteriorLightsElectricity(xi:xf);
% eagers = B_loads.InteriorLighting;
% [errors(figNum),Rsquare(figNum)] = elec_load_plot_setup(figTitle, figNum, Xe, ePlus, DofY, eagers, dt, dtE,ylab);
% 
% %compare exterior lighting
% figTitle = 'Exterior Lighting';
% figNum = 3;
% ylab = 'Electric Load (kW)';
% ePlus = mtr.ExteriorLightsElectricity(xi:xf);
% eagers = B_loads.ExteriorLighting;
% [errors(figNum),Rsquare(figNum)] = elec_load_plot_setup(figTitle, figNum, Xe, ePlus, DofY, eagers, dt, dtE,ylab);

%compare cooling electric loads
figTitle = 'Cooling';
figNum = 4;
ylab = 'Cooling (thermal kW)';
ePlus = mtr.CoolingElectricity(xi:xf);
eagers = Cooling / building.VariableStruct.COP_C;
[errors(figNum),Rsquare(figNum)] = elec_load_plot_setup(figTitle, figNum, Xe, ePlus, DofY, eagers, dt, dtE,ylab);

%compare heating electric loads
figTitle = 'Heating';
figNum = 5;
ylab = 'Heating (thermal kW)';
ePlus = mtr.HeatingElectricity(xi:xf)+mtr.HeatingGas(xi:xf);
eagers = Heating / building.VariableStruct.COP_H;
[errors(figNum),Rsquare(figNum)] = elec_load_plot_setup(figTitle, figNum, Xe, ePlus, DofY, eagers, dt, dtE,ylab);

%compare fan loads
figTitle = 'Fans';
figNum = 6;
ylab = 'Electric Load (kW)';
ePlus = mtr.FansElectricity(xi:xf);
eagers = Fan_Power;
[errors(figNum),Rsquare(figNum)] = elec_load_plot_setup(figTitle, figNum, Xe, ePlus, DofY, eagers, dt, dtE,ylab);

%compare net HVAC loads
figTitle = 'Net HVAC';
figNum = 7;
ylab = 'Energy (kW)';
ePlus = mtr.ElectricityHVAC(xi:xf)+mtr.HeatingGas(xi:xf);
eagers = B_loads.HVAC_electric;
[errors(figNum),Rsquare(figNum)] = elec_load_plot_setup(figTitle, figNum, Xe, ePlus, DofY, eagers, dt, dtE,ylab);

%compare net Building electric loads
figTitle = 'Net Building';
figNum = 8;
ylab = 'Electric Load (kW)';
ePlus = mtr.ElectricityNetFacility(xi:xf);
eagers = Total;
[errors(figNum),Rsquare(figNum)] = elec_load_plot_setup(figTitle, figNum, Xe, ePlus, DofY, eagers, dt, dtE,ylab);
% 
% %compare zone temp
% load('Tbuild.mat');
% figTitle = 'Zone Temperature';
% ylab = 'Temperature (C)';
% figNum = 9;
% ePlus = Tbuild;
% eagers = Tzone(2:end);
% % eagers = T_profile(2:end,1);
% [errors(figNum),Rsquare(figNum)] = elec_load_plot_setup(figTitle, figNum, DofY, ePlus, DofY, eagers, dt, dt,ylab);
% plotStep(DofY,Twall(2:end),'g')
% legend('EnergyPlus','Tzone','Twall')
% 
%     %compare air temp
%     figTitle = 'Zone Temperature';
%     fig = figure(10);
%     set(fig,'name',figTitle)
%     hold off
%     plotStep(date,Tzone(2:end),'b')
%     hold on
%     plotStep(date,Tzone_S(2:end),'r')
%     xlabel('Day of Year')
%     ylabel('Temperature (C)')
%     legend({'Forecast','Simulate'})
%     errors(1) = 100*(sum(Tzone_S.*dt)-sum(Tzone.*dt))/sum(Tzone.*dt);
%     fprintf('%s percent error:\t%f\n', figTitle, errors(1))
%     SSE = zeros(length(dt),1);
%     for i = 1:1:length(dt)
%         SSE(i) = (sum(Tzone_S((i-1)+1:i).*dt((i-1)+1:i)) - Tzone(i).*dt(i))^2;
%     end
%     Rsquare(1) = 1 - sum(SSE)/sum((Tzone.*dt-mean(Tzone.*dt)).^2);
%     fprintf('%s R^2 error:\t%f\n', figTitle, Rsquare(1))

end%Ends function CompareForecastSimulation

function plotStep(varargin)
Xe = varargin{1};
Ye = varargin{2};
if length(varargin) > 2
    lineSpec = {varargin{3:end}};
else
    lineSpec = {};
end
dt = Xe(2:end) - Xe(1:end-1);
[X,I] = sort([Xe(1)-dt(1);Xe;Xe(2:end)-dt+1e-9]);
Y = [Ye(1);Ye;Ye(2:end)];
Y = Y(I);
plot(X, Y, lineSpec{:})
end % Ends function plotStep

function [PercError,Rsquare] = elec_load_plot_setup(figTitle, figNum, Xe, E_Plus,DofY, Eagers, dt, dtE,ylab)
fig = figure(figNum);
set(fig,'name',figTitle)
hold off
plotStep(Xe,E_Plus,'b')
hold on
plotStep(DofY,Eagers,'r')
xlabel('Day of Year')
ylabel(ylab)
legend({'EnergyPlus','EAGERS'})
PercError = 100*(sum(Eagers.*dt)-sum(E_Plus.*dtE))/sum(E_Plus.*dtE);
fprintf('%s percent error:\t%f\n', figTitle, PercError)
n = round(dtE(1)/dt(1));
SSE = zeros(length(dtE),1);
for i = 1:1:length(dtE)
    SSE(i) = (sum(Eagers((i-1)*n+1:i*n).*dt((i-1)*n+1:i*n)) - E_Plus(i).*dtE(i))^2;
end
Rsquare = 1 - sum(SSE)/sum((E_Plus.*dtE-mean(E_Plus.*dtE)).^2);
fprintf('%s R^2 error:\t%f\n', figTitle, Rsquare)
end % Ends function elec_load_plot_setup

