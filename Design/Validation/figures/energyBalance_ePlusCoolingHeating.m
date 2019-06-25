%ENERGYBALANCE_EPLUSCOOLINGHEATING

%% Constants
BUILDING = 'MediumOffice';

%% Load EPlus zone temps, cooling, heating, date
thisFilePath = mfilename('fullpath');
splitted = strsplit(thisFilePath, filesep);
thisFileDir = strjoin({splitted{1:end-1}}, filesep);
mtrFile = fullfile(pwd, '..', strcat(BUILDING,'.mat'));
mztFile = fullfile(pwd, '..', strcat('meanZoneTemps_',BUILDING,'.mat'));
load(mtrFile);
load(mztFile);
cooling = mtr.CoolingElectricity;
heating = mtr.HeatingElectricity;
date = mtr.Timestamp;

%% Energy balance
Tact = Tact + ((Tdb(t) - Tact)*Build.Area/Build.VariableStruct.Resistance + ...
    InternalGains(t) + (Tsupply - Tact)*Cp_Air*Flow)*dt(t)/ ...
    (Build.VariableStruct.Capacitance*Build.Area); %net change in temperature

%% Clean up extraneous variables
clear mtrFile mztFile
