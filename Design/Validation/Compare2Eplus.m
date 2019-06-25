function [errors,Rsquare] = Compare2Eplus(building,weather,Date)
% Compare a forecasted building to it's energyPlus counterpart
% building is the EAGERS building structure, weather is a structure of dry
% bulb, wetbuld and relative humididty, and solar irradiation,
% Date is the timsteps to be compared (must be sequential).

%% Handle input DATE
[s1, ~] = size(Date);
if s1 == 1
    Date = Date';
end
Location = struct('Latitude',40, 'Longitude',-105, 'TimeZone',-7);
%% Load EPlus meter file
global Model_dir
load(fullfile(Model_dir,'Design','Validation',strcat(building.Name,'.mat')))% building
%% Run EAGERS method
Tzone = 17.4;
Twall = 17;
iWeather = interpolate_weather(weather,Date);
Tdb = iWeather.Tdb;
RH = iWeather.RH;
SG = solar_gain(building,Date,Location,iWeather);
ExternalGains = SG.Walls + SG.Roof;
B_loads = building_loads(building,Date,SG);
xF = nnz(Date<Date(1)+7);
wuDate = Date(1:xF);
wuWeather = interpolate_weather(weather,wuDate);
wuSG = solar_gain(building,wuDate,Location,wuWeather);
wuLoads = building_loads(building,wuDate,wuSG);
wuExternalGains = wuSG.Walls + wuSG.Roof;
[~,~,~,Tzone,Twall,~] = building_profile(building,wuDate,wuLoads.InternalGains,wuExternalGains,wuWeather.Tdb,wuWeather.RH,Tzone(end,1),Twall(end,1));

% [Cooling, Heating, Fan_Power,Damper,T_profile] = building_profile2(building,Date,B_loads.InternalGains,ExternalGains,Tdb,RH,[17.5,15,11,5]);
[Cooling, Heating, Fan_Power,Tzone,Twall,Damper] = building_profile(building,Date,B_loads.InternalGains,ExternalGains,Tdb,RH,Tzone(end,1),Twall(end,1));
B_loads.HVAC_electric = Heating/building.VariableStruct.COP_H + Cooling/building.VariableStruct.COP_C + Fan_Power;
Total = B_loads.Equipment + B_loads.InteriorLighting + B_loads.ExteriorLighting + B_loads.HVAC_electric +  B_loads.OtherLoads;

D = datevec(Date(1));
DofY = Date - datenum([D(1), 1,1]);
dt = Date(2:end) - Date(1:end-1);
dt = [dt(1);dt];
if isempty(mtr.Timestamp)
    EplusTime = linspace(datenum([2017,1,1,1,0,0]),datenum([2018,1,1]),8760)';
else
    EplusTime = datenum(datevec(mtr.Timestamp));
end
D = datevec(EplusTime(1));
EplusTime = EplusTime - datenum([D(1), 1,1]);

xi = max(1,nnz((EplusTime)<=DofY(1)));
xf = nnz((EplusTime)<=DofY(end));
Xe = EplusTime(xi:xf);
if xi == 1
    dtE = Xe - ([0;EplusTime(1:xf-1)]);
else
    dtE = Xe - (EplusTime(xi-1:xf-1));
end

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

%compare zone temp
load('Tbuild.mat');
figTitle = 'Zone Temperature';
ylab = 'Temperature (C)';
figNum = 9;
ePlus = Tbuild;
eagers = Tzone(2:end);
% eagers = T_profile(2:end,1);
[errors(figNum),Rsquare(figNum)] = elec_load_plot_setup(figTitle, figNum, DofY, ePlus, DofY, eagers, dt, dt,ylab);
plotStep(DofY,Twall(2:end),'g')
legend('EnergyPlus','Tzone','Twall')
% plotStep(DofY,T_profile(2:end,2),'g')
% plotStep(DofY,T_profile(2:end,3),'c')
% figure
% plot(Date,Tbuild,'g')
% hold on
% plot(Date,Tzone,'k')
% plot(Date,Twall,'c')


% 
% figure
% plotStep(Xe(351*24+1:362*24),mtr.HeatingElectricity(351*24+1:362*24)*3.7382,'b')
% hold on
% plotStep(DofY(351*24+1:362*24),Heating(351*24+1:362*24),'r')
% xlabel('Day of Year')
% ylabel('Heating Load (kW)')
% legend({'EnergyPlus','EAGERS'})
% 
% figure
% plotStep(Xe(169*24+1:181*24),mtr.CoolingElectricity(169*24+1:181*24)*3.6504,'b')
% hold on
% plotStep(DofY(169*24+1:181*24),Cooling(169*24+1:181*24),'r')
% xlabel('Day of Year')
% ylabel('Cooling Load (kW)')
% legend({'EnergyPlus','EAGERS'})
% 
% figure
% plotStep(Xe(351*24+1:362*24),mtr.ElectricityNetFacility(351*24+1:362*24)-mtr.ElectricityHVAC(351*24+1:362*24),'b')
% hold on
% plotStep(DofY(351*24+1:362*24),Total(351*24+1:362*24)-HVAC_electric(351*24+1:362*24),'r')
% xlabel('Day of Year')
% ylabel('Uncontrolled Electric Load (kW)')
% legend({'EnergyPlus','EAGERS'})
end % Ends function Compare2Eplus


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
