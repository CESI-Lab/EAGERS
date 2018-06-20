function sg = solar_gain(building,date,location,weather)
%% Solar gain
n_s = length(date);
or = building.VariableStruct.Orientation;
sg_diffuse = building.VariableStruct.WallArea*weather.irradDiffHorz;%Diffuse irradiance (W)
[sg.Sunrise, sg.Sunset, azimuth, zenith] = solar_calc(location.Longitude, location.Latitude, location.TimeZone, date);
dir_norm = [cosd(azimuth).*cosd(90 - zenith), sind(azimuth).*cosd(90 - zenith),sind(90 - zenith)];%Direct normal (incoming vector of sunlight)
n = ones(n_s,1);
sg_direct = 0.25*building.VariableStruct.WallArea*weather.irradDireNorm.*(max(0,dot(n*[cosd(or),sind(or),0],dir_norm,2)) + max(0,dot(n*[-cosd(or),-sind(or),0],dir_norm,2)) + max(0,dot(n*[-sind(or),cosd(or),0],dir_norm,2)) + max(0,dot(n*[sind(or),-cosd(or),0],dir_norm,2)));
sg.Windows = building.VariableStruct.WindowTransmittance*building.VariableStruct.WindowWallRatio*(sg_direct+sg_diffuse)/1000;%solar gain (heat) through windows http://www.commercialwindows.org/vt.php
sg.Walls = building.VariableStruct.WallAbsorption*(1-building.VariableStruct.WindowWallRatio)*(sg_direct+sg_diffuse)/1000;%solar gain (heat) absorbed by walls
sg.VisibleLight = building.VariableStruct.LightTransmittance*building.VariableStruct.WindowWallRatio*(sg_direct+sg_diffuse)/1000;%visible light transmitted through windows
sg.Roof = building.VariableStruct.RoofArea*building.VariableStruct.WallAbsorption*weather.irradDireNorm.*max(0,dot(n*[0,0,1],dir_norm,2))/1000;%solar gain (heat) through roof
end%Ends function solar_gain