function wind_speed = wind_speed_calc(wind,height,terrain)
d_met = 270;%boundary layer thickness at meterological site
z_met = 10;%elevation of met site
a_met = 0.14; %coefficient at met site
switch terrain
    case 'City'
        a = .33;
        d = 460;
    case 'Urban'
        a = .22;
        d = 370;
    case 'Ocean'
        a = 0.1;
        d = 210;
    case 'Flat'
        a = 0.14;
        d = 270;
    case 'Rough'
        a = 0.22;
        d = 370;
    otherwise
        disp('what is terrain')
end
wind_speed = wind*(d_met/z_met)^a_met*(height/d).^a;
end%Ends function wind_speed_calc