Model_dir = strrep(which('building_validate.m'),'building_validate.m','');
addpath(fullfile(Model_dir,'Optimization','BasicFunctions'));
addpath(fullfile(Model_dir,'Optimization','BuildingModeling'));
addpath(fullfile(Model_dir,'System Library','Weather'));
addpath(fullfile(Model_dir,'System Library','Buildings','E_PLUS'));

weather = import_weather('4A_USA_MD_BALTIMORE_TMY2.epw');


b_names = {'RefBldgFullServiceRestaurantNew2004_v1';...%1
            'RefBldgHospitalNew2004_v1';...%2
            'RefBldgLargeHotelNew2004_v1';...%3
            'RefBldgLargeOfficeNew2004_v1';...%4
            'RefBldgMediumOfficeNew2004_v1';...%5
            'RefBldgMidriseApartmentNew2004_v1';...%6
            'RefBldgOutPatientNew2004_v1';...%7
            'RefBldgPrimarySchoolNew2004_v1';...%8
            'RefBldgQuickServiceRestaurantNew2004_v1';...%9
            'RefBldgSecondarySchoolNew2004_v1';...%10
            'RefBldgSmallHotelNew2004_v1';...%11
            'RefBldgSmallOfficeNew2004_v1';...%12
            'RefBldgStand-aloneRetailNew2004_v1';...%13
            'RefBldgStripMallNew2004_v1';...%14
            'RefBldgSuperMarketNew2004_v1';...%15
            'RefBldgWarehouseNew2004_v1';};%16

name = b_names{1};
building = import_idf(strcat(name,'.4_7.2_4A_USA_MD_BALTIMORE'));
size_info = import_sizes(strcat(name,'.4_7.2_4A_USA_MD_BALTIMORE'),building.zones.name,building.HVAC.loop.name,building.zones.multiplier);
building = size_hvac(building,weather,size_info);
date = building.sim_date;
% index = 1:240;
% index = 5185:5424;
index = 1:8760;
date2 = date([index,index(end)+1]);%shorten simulation
[T_zone,T_surf,m_v_zone,m_sys,Q_transfer,heat_elec,heat_gas,cooling,Q_test_plots,T_windows_int] = ...
    zone_simulate(building,[],[],weather,date2,[]);

profile = load_sched(building.schedule,building.holidays,date2(2:end),[],[]);
weather2 = interpolate_weather(weather,date2(2:end));
dt = (date2(2:end)-date2(1:end-1))*24*3600;
T_air_z = weather2.DrybulbC*ones(1,length(building.ctf.z_height)) - ones(length(index),1)*0.0065*building.ctf.z_height';
P = weather2.PressurePa/1000; % atmospheric pressure (kPa)
T_dp = 273 + weather2.DewpointC;
P_H2O_dp = exp((-5.8002206e3)./T_dp + 1.3914993 - 4.8640239e-2*T_dp + 4.1764768e-5*T_dp.^2 - 1.4452093e-8*T_dp.^3 + 6.5459673*log(T_dp))/1000; %saturated water vapor pressure at dewpoint
m_v_air = .621945*(P_H2O_dp./(P-P_H2O_dp));%mass fraction of water in air 
frost = zeros(length(date2),length(building.cases.name)); 
[net_building,loads,occupancy,mixing,infiltration,Q_zone_gain,Q_surf_absorb,...
    W_zone_gain,internal_irrad_window,T_set,cos_phi_windows,frost(2:end,:)] = ...
    zone_loads(building,date2(2:end),profile,weather2.DNIWm2,weather2.DHIWm2,T_zone(2:end,:),m_v_zone(2:end,:),dt,T_air_z,m_v_air,frost(1:end-1,:));

%% import EnergyPlus results as table
e_plus = readtable(strcat(name,'.4_7.2_4A_USA_MD_BALTIMORE.csv'));
e_plus.(e_plus.Properties.VariableNames{end}) = str2double(e_plus.(e_plus.Properties.VariableNames{end}));
zonenames = building.zones.name;
lc = {'a';'b';'c';'d';'e';'f';'g';'h';'i';'j';'k';'l';'m';'n';'o';'p';'q';'r';'s';'t';'u';'v';'w';'x';'y';'z';};
uc = {'A';'B';'C';'D';'E';'F';'G';'H';'I';'J';'K';'L';'M';'N';'O';'P';'Q';'R';'S';'T';'U';'V';'W';'X';'Y';'Z';};
for i = 1:1:length(zonenames)
    for j = 1:1:26
        zonenames(i) = {strrep(zonenames{i},lc{j},uc{j})};
    end
%     zonenames(i) = {strrep(zonenames{i},'_','')};
end
plot_z = [1;2;3;];

d1 = date(1);
figure(1)
plot(date2(2:end)-d1,heat_gas/1000,'k');
hold on
plot(date2(2:end)-d1,e_plus.Heating_Gas_J__Hourly_(index)/3600/1000,'r');
legend({'WSU';'EnergyPlus';})
xlabel('Day of year')
ylabel('heating (kW)')
figure(2)
plot(date2(2:end)-d1,cooling/1000,'k');
hold on
plot(date2(2:end)-d1,e_plus.Cooling_Electricity_J__Hourly_(index)/3600/1000,'r');
legend({'WSU';'EnergyPlus';})
xlabel('Day of year')
ylabel('cooling (kW)')

figure(3)
plot(date2(2:end)-d1,m_v_zone(2:end,:));
legend(building.zones.name)
xlabel('Day of year')
ylabel('zone moisture content (kg H2O/kg air)')

figure(4)
hold on
T_interior = T_surf(:,building.ctf.sur_state1);
for p = 1:1:length(plot_z)
    s_area = building.surfaces.area(building.zones.surfaces{plot_z(p)});
    T_interior_avg = T_interior(2:end,building.zones.surfaces{plot_z(p)})*s_area/sum(s_area);
    plot(date2(2:end)-d1,T_interior_avg);
end
legend(building.zones.name(plot_z))
xlabel('Day of year')
ylabel('zone surface temperature (C)')

p_n = 4;
for p = 1:1:length(plot_z)
    figure(p_n+p)
    plot(date2(2:end)-d1,max(0,-Q_transfer(:,plot_z(p)))/1000,'k');
    hold on
    field = strcat('Cooling_EnergyTransfer_Zone_',zonenames{plot_z(p)},'_J__Hourly_');
    field = field(1:min(63,length(field)));
    plot(date2(2:end)-d1,e_plus.(field)(index)/3600/1000,'r');
    legend({'WSU';'EnergyPlus';})
    xlabel('Day of year')
    ylabel(strcat(zonenames{plot_z(p)},' zone cooling (kW)'))
end

p_n = p_n + length(plot_z);
for p = 1:1:length(plot_z)
    figure(p_n+p)
    plot(date2(2:end)-d1,T_zone(2:end,plot_z(p)),'k');
    hold on
    field = strcat(zonenames{plot_z(p)},'_ZoneAirTemperature_C__Hourly_');
    field = field(1:min(63,length(field)));
    plot(date2(2:end)-d1,e_plus.(field)(index),'r');
    legend({'WSU';'EnergyPlus';})
    xlabel('Day of year')
    ylabel(strcat(zonenames{plot_z(p)},' zone temperature (C)'))
end

% %don't plot all internal gains so choose from list
% g_list = {'LatentGain';'RadiantHeating';'VisibleRadiationHeating';'ConvectiveHeating';};
% g_list2 = {'latent';'radiant';'visible';'convected';};
% k = 4;
% p_n = p_n + length(plot_z);
% for p = 1:1:length(plot_z)
%     figure(p_n+p)
%     plot(date2(2:end)-d1,loads.(g_list2{k})(:,plot_z(p))/1000,'k');
%     hold on
%     field = strcat(zonenames{plot_z(p)},'_ZoneTotalInternal',g_list{k},'Rate_W__Hourly_');
%     field = field(1:min(63,length(field)));
%     plot(date2(2:end)-d1,e_plus.(field)(index)/1000,'r');
%     legend({'WSU';'EnergyPlus';})
%     xlabel('Day of year')
%     ylabel(strcat(zonenames{plot_z(p)},' zone ',g_list2{k},' gain (kW)'))
% end
% 
% p_n = p_n + length(plot_z);
% for p = 1:1:length(plot_z)
%     figure(p_n+p)
%     plot(date2(2:end)-d1,Q_test_plots.Infiltration(:,plot_z(p))/1000,'k');
%     hold on
%     field = strcat(zonenames{plot_z(p)},'_ZoneAirHeatBalanceOutdoorAirTransferRate_W__Hourly_');
%     field = field(1:min(63,length(field)));
%     plot(date2(2:end)-d1,e_plus.(field)(index)/1000,'r');
%     legend({'WSU';'EnergyPlus';})
%     xlabel('Day of year')
%     ylabel(strcat(zonenames{plot_z(p)},' zone outside air gain (kW)'))
% end
% 
% p_n = p_n + length(plot_z);
% for p = 1:1:length(plot_z)
%     figure(p_n+p)
%     plot(date2(2:end)-d1,Q_test_plots.surface_conv(:,plot_z(p))/1000,'k');
%     hold on
%     field = strcat(zonenames{plot_z(p)},'_ZoneAirHeatBalanceSurfaceConvectionRate_W__Hourly_');
%     field = field(1:min(63,length(field)));
%     plot(date2(2:end)-d1,e_plus.(field)(index)/1000,'r');
%     legend({'WSU';'EnergyPlus';})
%     xlabel('Day of year')
%     ylabel(strcat(zonenames{plot_z(p)},' zone surface convection gain (kW)'))
% end
% 
% p_n = p_n + length(plot_z);
% for p = 1:1:length(plot_z)
%     figure(p_n+p)
%     plot(date2(2:end)-d1,Q_test_plots.air_balance_convective(:,plot_z(p))/1000,'k');
%     hold on
%     field = strcat(zonenames{plot_z(p)},'_ZoneAirHeatBalanceInternalConvectiveHeatGainRate_W__Hour');
%     field = field(1:min(63,length(field)));
%     plot(date2(2:end)-d1,e_plus.(field)(index)/1000,'r');
%     legend({'WSU';'EnergyPlus';})
%     xlabel('Day of year')
%     ylabel(strcat(zonenames{plot_z(p)},' zone air balance convection gain (kW)'))
% end