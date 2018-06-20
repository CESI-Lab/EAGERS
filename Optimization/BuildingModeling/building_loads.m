function B_Loads = building_loads(Build,Date,SolarGain)
% Time
[Y,Month,D,H,M,S] = datevec(Date);
h1 = H(1) + M(1) + S(1);
D(H==0&M==0&S==0) = D(H==0&M==0&S==0) - 1;
if h1 == 0 % put start date back.
    D(1) = D(1) + 1;
end
days = ceil(Date(end) - datenum([Y(1),Month(1),D(1)]));
sd = floor(Date(1)); % start date
nS = length(Date); % number timesteps

%% Fill out schedule values for each timestep in input Date
% If next point in schedule is shorter than ramp, then ramp shortens to half the gap
p = 1;
DateIndex = cell(days);
for d = 1:1:days
    p2 = p;
    while p2<nS && Date(p2+1)<=sd+d
        p2 = p2+1;
    end
    DateIndex(d) = {p:p2};
    p = p2+1;
end
% Schedules
sched = fieldnames(Build.Schedule);
for i = 1:1:length(sched)
    if ~isfield(Build.Schedule.(sched{i}),'Timestamp') %use EnergyPlus stype weekday, Sat, sun schedules
        Profile.(sched{i}) = load_sched(Build,Date,sched{i}); 
    else %use data saved in schedule as representative random schedule behavior (at least 1 day of data is necessary)
        Profile.(sched{i}) = load_sched2(Build.Schedule.(sched{i}),Date); 
    end
end
%% Equipment load values
B_Loads.Equipment = Build.Area*Build.VariableStruct.equipment/1000*Profile.equipment; % kW of equipment load

%% Other load values
if isfield(Build.VariableStruct,'DataCenter')%data center is an equipment load but not an internal gain into the zones
    B_Loads.Equipment = B_Loads.Equipment + Build.VariableStruct.DataCenter*Profile.datacenter;
end
B_Loads.OtherLoads = zeros(nS,1);
if isfield(Build.VariableStruct,'WaterSystem')
    B_Loads.OtherLoads = Build.Area*Build.VariableStruct.WaterSystem*Profile.watersystem;
end
B_Loads.DCloads = [];
if isfield(Build.VariableStruct,'DCloads')
    B_Loads.DCloads = Build.Area*Build.VariableStruct.DCloads*Profile.DCloads;
end

%% Exterior lighting load values
if isfield(Profile,'exteriorlights_solarcontroled')
    ExtLights_B = .04;%.125; % fraction of day exterior lights remain on after sunrise
    for day = 1:days
        currDate = floor(Date(DateIndex{day}(1)));
        HoD = 24*(Date(DateIndex{day}) - currDate); % hour of day
        frac_dark = ones(length(HoD),1);
        Sunrise = SolarGain.Sunrise(DateIndex{day}(1));
        Sunset = SolarGain.Sunset(DateIndex{day}(1));

        % Find hour at which sunrise occurs
        h_sr = 1;
        while HoD(h_sr)<(Sunrise*24+ExtLights_B) && h_sr<length(HoD)
            h_sr = h_sr+1;
        end
        if h_sr == 1
            frac_dark(1) = (Sunrise*24+ExtLights_B)/HoD(1);
        else
            frac_dark(h_sr) = ((Sunrise*24+ExtLights_B) - HoD(h_sr-1))/(HoD(h_sr)-HoD(h_sr-1));
        end
        
        % Find hour at which sunset occurs
        h_ss = h_sr+1;
        if h_ss < length(HoD)
            while HoD(h_ss)<(Sunset*24-ExtLights_B)  && h_ss<length(HoD)
                h_ss = h_ss+1;
            end
            frac_dark(h_ss) = 1-((Sunset*24-ExtLights_B) - HoD(h_ss-1))/(HoD(h_ss)-HoD(h_ss-1));
        end
        frac_dark(h_sr+1:h_ss-1) = 0;
        Profile.exteriorlights_solarcontroled(DateIndex{day}) = min(Profile.exteriorlights_solarcontroled(DateIndex{day}'),frac_dark);
    end
    B_Loads.ExteriorLighting = Build.VariableStruct.ExteriorLights*Profile.exteriorlights_fixed + Build.VariableStruct.ExteriorLights*Profile.exteriorlights_solarcontroled;
else
    B_Loads.ExteriorLighting = Build.VariableStruct.ExteriorLights*Profile.exteriorlights;
end

%% Interior lighting load values
B_Loads.InteriorLighting = Build.Area*Build.VariableStruct.InteriorLights/1000*Profile.interiorlights;% kW of lighting load
if any(SolarGain.VisibleLight>0)
    B_Loads.InteriorLighting = B_Loads.InteriorLighting.*(1 - Build.VariableStruct.DaylightPercent*SolarGain.VisibleLight/max(SolarGain.VisibleLight)); %daylighting
end

%% Internal gains
internalGainOccupants = Build.Area * Build.VariableStruct.occupancy * Profile.occupancy * 0.120; % heat from occupants (120 W)
B_Loads.InternalGains = internalGainOccupants + B_Loads.Equipment + B_Loads.InteriorLighting + SolarGain.Windows;
end%Ends function building_loads
