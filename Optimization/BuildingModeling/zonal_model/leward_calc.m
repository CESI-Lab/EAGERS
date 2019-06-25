function leward = leward_calc(wind,surf_normal)
threshold = 100; %degrees difference between wind and surface normal to be considered leward
surf = atand(surf_normal(:,1)./surf_normal(:,2));
quad2 = surf_normal(:,1)>=0 & surf_normal(:,2)<=0;
quad3 = surf_normal(:,1)<0 & surf_normal(:,2)<=0;
surf(quad2 | quad3) = surf(quad2 | quad3) + 180;
surf(surf_normal(:,2)==0 & quad2) = 90;
surf(surf_normal(:,2)==0 & quad3) = 270;
wind2 = wind-360;
leward = abs(wind-surf)>threshold & abs(wind2-surf)>threshold;
% leward = (cosd(weather.Wdirdegrees)*windows.normal(:,2) + sind(weather.Wdirdegrees)*windows.normal(:,1))<0;%works for threshold of 90 degrees
end%Ends function leward_calc
