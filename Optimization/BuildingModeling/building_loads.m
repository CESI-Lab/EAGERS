function loads = building_loads(building,date,solar_gain)
n_s = length(date); % number timesteps
%% Fill out schedule values for each timestep in input Date
% If next point in schedule is shorter than ramp, then ramp shortens to half the gap
% Schedules
profile = load_sched(building.Schedule,date,[]); 
%% Equipment load values
loads.Equipment = building.Area*building.VariableStruct.equipment/1000*profile.equipment; % kW of equipment load

%% Other load values
if isfield(building.VariableStruct,'DataCenter')%data center is an equipment load but not an internal gain into the zones
    loads.Equipment = loads.Equipment + building.VariableStruct.DataCenter*profile.datacenter;
end
loads.OtherLoads = zeros(n_s,1);
if isfield(building.VariableStruct,'WaterSystem')
    loads.OtherLoads = building.Area*building.VariableStruct.WaterSystem*profile.watersystem;
end
loads.DCloads = [];
if isfield(building.VariableStruct,'DCloads')
    loads.DCloads = building.Area*building.VariableStruct.DCloads*profile.DCloads;
end

%% Exterior lighting load values
if isfield(profile,'exteriorlights_solarcontroled')
    loads.ExteriorLighting  = astronomical_clock_controlled_exterior_lights(date,solar_gain,building.VariableStruct.ExteriorLights,profile.exteriorlights_solarcontroled,profile.exteriorlights_fixed);
else
    loads.ExteriorLighting = building.VariableStruct.ExteriorLights*profile.exteriorlights;
end

%% Interior lighting load values
loads.InteriorLighting = building.Area*building.VariableStruct.InteriorLights/1000*profile.interiorlights;% kW of lighting load
if any(solar_gain.VisibleLight>0)
    loads.InteriorLighting = loads.InteriorLighting.*(1 - building.VariableStruct.DaylightPercent*solar_gain.VisibleLight/max(solar_gain.VisibleLight)); %daylighting
end

%% Internal gains
internalGainOccupants = building.Area * building.VariableStruct.occupancy * profile.occupancy * 0.120; % heat from occupants (120 W)
loads.InternalGains = internalGainOccupants + loads.Equipment + loads.InteriorLighting + solar_gain.Windows;
end%Ends function building_loads

function lights = astronomical_clock_controlled_exterior_lights(date,solar_gain,max_light,sched1,sched2)
n_s = length(date);
shift = .04;%.125; % fraction of day exterior lights remain on after sunrise
frac_dark = ones(n_s,1);
if n_s == 1
    dt = [1,1];
else
    dt = 24*(date(2:end) - date(1:end-1));
    dt = [dt(1);dt];%shift things by 1 so frac_dark at first date can be calculated
end
date = [date(1)-dt(1);date];
hour =  24*(date - floor(date)); % hour of day
for t = 2:1:n_s+1
    if hour(t)<=(solar_gain.Sunrise(t-1)*24+shift)
        %still dark
    elseif hour(t)>(solar_gain.Sunrise(t-1)*24+shift) && hour(t-1)<=(solar_gain.Sunrise(t-1)*24+shift)%sunrise
        frac_dark(t-1) = 1-((solar_gain.Sunrise(t-1)*24+shift) - hour(t-1))/dt(t-1);
    elseif hour(t)<(solar_gain.Sunset(t-1)*24-shift)
        frac_dark(t-1) = 0;
    elseif hour(t)>(solar_gain.Sunset(t-1)*24-shift) && hour(t-1)<=(solar_gain.Sunset(t-1)*24-shift)%sunset
        frac_dark(t-1) = ((solar_gain.Sunset(t-1)*24-shift) - hour(t-1))/dt(t-1);
    end
end
lights = max_light*sched2 + max_light*min(sched1,frac_dark);
end%Ends function astronomical_clock_controlled_exterior_lights
