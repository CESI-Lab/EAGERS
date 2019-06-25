%OVERLAYDEMANDTEMPS_MEANZONE_ANDMEANWTDBLDG
%Make 2 plots:
%1.     Cooling and heating overlay mean zone temperatures
%2.     Cooling and heating overlay area-weighted mean total building
%       temperature

%% Constants
BUILDING = 'MediumOffice';

%% Load EPlus zone temps, cooling, heating, date
thisFilePath = mfilename('fullpath');
splitted = strsplit(thisFilePath, filesep);
thisFileDir = strjoin({splitted{1:end-1}}, filesep);
mtrFile = fullfile(thisFileDir, '..', strcat(BUILDING,'.mat'));
mztFile = fullfile(thisFileDir, '..', strcat('meanZoneTemps_',BUILDING,'.mat'));
load(mtrFile);
load(mztFile);
cooling = mtr.CoolingElectricity;
heating = mtr.HeatingElectricity;
date = mtr.Timestamp;

%% Plot mean zone temperatures
figure
stairs(date, cooling)
hold on
stairs(date, heating)
fdNames = fieldnames(meanZoneTemps);
for i = 1:1:length(fdNames)
    stairs(date, meanZoneTemps.(fdNames{i}))
end
hold off
set(gcf, 'Position', [408 1.096e+03 1144 420])
legend({'Cooling','Heating','Zone Temperatures'})
title({'EnergyPlus Cooling and Heating Demand','overlay Mean Zone Temperatures'})
xlabel('Date')
ylabel({'Demand [kW]','Temperature [C]'})

%% Plot area-weighted mean total building temperature
figure
stairs(date, cooling)
hold on
stairs(date, heating)
stairs(date, meanZoneTemps.BUILDING_WEIGHTED_AVG, 'k')
hold off
set(gcf, 'Position', [408, 1.096e+03, 1144, 420])
legend({'Cooling Demand','Heating Demand','Area-Wtd TotalMean Bldg Temp'})
title({'EnergyPlus Cooling and Heating Demand', ...
    'overlay Area-Weighted Mean Total Building Temperatures'})
xlabel('Date')
ylabel({'Demand [kW]','Temperature [C]'})

%% Clean up extraneous variables
clear mtrFile mztFile fdNames i
