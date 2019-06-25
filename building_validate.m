function [building,loads,gains,model_fit,facility,e_plus] = building_validate(b_num,plots_on)
%% [building,loads,gains,model_fit,facility,e_plus] = building_validate(1:16,0);
%% [building,loads,gains,model_fit,facility,e_plus] = building_validate(4,1);
Model_dir = strrep(which('building_validate.m'),'building_validate.m','');
addpath(Model_dir);
addpath(fullfile(Model_dir,'Optimization','BasicFunctions'));
addpath(fullfile(Model_dir,'Optimization','BuildingModeling','zonal_model'));
addpath(fullfile(Model_dir,'System_Library','Weather'));
addpath(fullfile(Model_dir,'System_Library','Buildings','E_PLUS'));

weather = import_weather('4A_USA_MD_BALTIMORE_TMY2.epw');

b_names = {'FullServiceRestaurant';...%1
            'Hospital';...%2
            'LargeHotel';...%3
            'LargeOffice';...%4
            'MediumOffice';...%5
            'MidriseApartment';...%6
            'OutPatient';...%7
            'PrimarySchool';...%8
            'QuickServiceRestaurant';...%9
            'SecondarySchool';...%10
            'SmallHotel';...%11
            'SmallOffice';...%12
            'Stand-aloneRetail';...%13
            'StripMall';...%14
            'SuperMarket';...%15
            'Warehouse';};%16
        
EP = {'Heating_Gas_J__Hourly_';'Heating_Electricity_J__Hourly_';'WaterSystems_Gas_J__Hourly_';'Cooling_Electricity_J__Hourly_';'HeatRejection_Electricity_J__Hourly_';'Fans_Electricity_J__Hourly_';'Pumps_Electricity_J__Hourly_';};
EAG = {'heat_gas';'heat_elec';'water_gas';'cool_elec';'tower_elec';'fan_elec';'pumps';};
label_x = {'heating --- gas (kW)';'heating --- electric (kW)';'water system --- gas (kW)';'cooling (kW)';'cooling tower fans (kW)';'HVAC fans (kW)';'water pumps --- (kW)';};

lc = {'a';'b';'c';'d';'e';'f';'g';'h';'i';'j';'k';'l';'m';'n';'o';'p';'q';'r';'s';'t';'u';'v';'w';'x';'y';'z';};
uc = {'A';'B';'C';'D';'E';'F';'G';'H';'I';'J';'K';'L';'M';'N';'O';'P';'Q';'R';'S';'T';'U';'V';'W';'X';'Y';'Z';};

for b = 1:1:16
    if any(b_num==b)
    tic
    building = import_idf(strcat('RefBldg',b_names{b},'New2004_v1.4_7.2_4A_USA_MD_BALTIMORE'));
    building = import_sizes(strcat('RefBldg',b_names{b},'New2004_v1.4_7.2_4A_USA_MD_BALTIMORE'),building);
    building.ctf = build_ctf(building.zones.name,building.surfaces,building.windows);
    
    date = building.sim_date;
    index = 1:240;
%     index = 2161:2640;
%     index = 3361:3840;
%     index = 5185:5424;
%     index = 1:8760;
    date2 = date([index,index(end)+1]);%shorten simulation
    [T_zone,T_surf,m_v_zone,m_sys,Q_transfer,facility,Q_test_plots] = ...
        run_zonal_building(building,[],[],weather,date2,[]);

    %% import EnergyPlus results as table
    e_plus = readtable(strcat('RefBldg',b_names{b},'New2004_v1.4_7.2_4A_USA_MD_BALTIMORE.csv'));
    zonenames = building.zones.name;
    for i = 1:1:length(zonenames)
        for j = 1:1:26
            zonenames(i) = {strrep(zonenames{i},lc{j},uc{j})};
        end
        zonenames(i) = {strrep(zonenames{i},' ','')};
        zonenames(i) = {strrep(zonenames{i},'-','_')};
    end
    d1 = date(1);
    F = fieldnames(e_plus);    
    for j = 1:1:length(EP)
        if any(strcmp(F,EP{j})) && any(e_plus.(EP{j})(index)>0)
            x = e_plus.(EP{j})(index)/3600/1000;
            y = facility.(EAG{j})/1000;
            COD = 1 - sum((y - x).^2)/sum((x-mean(x)).^2);%Coefficient of determination
            COD_B = 1 - sum((y - x).^2)/sum((y-mean(y)).^2);
            model_fit.(EAG{j}).Pearson(b) = sum((x-mean(x)).*(y - mean(y)))^2/(sum((x-mean(x)).^2)*sum((y-mean(y)).^2));%Pearson product moment correlation coefficient
            model_fit.(EAG{j}).COD(b) = max(COD,COD_B);
            model_fit.(EAG{j}).total_energy_percent(b) = sum(y)/sum(x);
        end
    end
    toc
    profile = load_sched(building.schedule,building.holidays,date2(2:end),[],[]);
    weather2 = interpolate_weather(weather,date2(2:end));
    dt = (date2(2:end)-date2(1:end-1))*24*3600;
    frost = zeros(length(date2),length(building.cases.name)); 
    [loads,gains,occupancy,mixing,infiltration,T_set,w_set,cos_phi_windows,frost(2:end,:),T_mains,water_heat] = ...
        zone_loads(building,date2(2:end),profile,weather2,T_zone(2:end,:),m_v_zone(2:end,:),dt,frost(1:end-1,:));
% 
%     net_building.interior_lighting = loads.lighting_internal*loads.multiplier;
%     net_building.equipment = loads.plug_load*loads.multiplier;
%     net_building.internal_gain = loads.internal_gain*loads.multiplier;
%     net_building.non_hvac_electric = net_building.interior_lighting + net_building.equipment + loads.exterior.lighting + loads.exterior.equipment + sum(loads.case.electric,2) + sum(loads.rack.electric,2);
    end
end
if length(b_num) ==1 && plots_on
    % plot_z = [1;2;3;];
    plot_z = (1:length(zonenames))';
    % plot_z = [2;3;4;5;];
    for j = 1:1:length(EP)
        if any(strcmp(F,EP{j})) && any(e_plus.(EP{j})(index)>0)
            x = e_plus.(EP{j})(index)/3600/1000;
            y = facility.(EAG{j})/1000;
            plot_fig(j,date2(2:end)-d1,x,y,label_x{j})
        end
    end

    % figure(length(EP)+1)
    % plot(date2(2:end)-d1,m_v_zone(2:end,:));
    % legend(building.zones.name)
    % xlabel('Day of year')
    % ylabel('zone moisture content (kg H2O/kg air)')
    % 
    % figure(length(EP)+2)
    % hold on
    % T_interior = T_surf(:,building.ctf.sur_state1);
    % for p = 1:1:length(plot_z)
    %     s_area = building.surfaces.area(building.zones.surfaces{plot_z(p)});
    %     T_interior_avg = T_interior(2:end,building.zones.surfaces{plot_z(p)})*s_area/sum(s_area);
    %     plot(date2(2:end)-d1,T_interior_avg);
    % end
    % legend(building.zones.name(plot_z))
    % xlabel('Day of year')
    % ylabel('zone surface temperature (C)')

    p_n = length(EP)+2;
    F = fieldnames(e_plus);
    for j = 1:1:length(F)-3
        if iscell(e_plus.(F{j}))
            e_plus.(F{j}) = str2double(e_plus.(F{j}));
        end
    end
    for p = 1:1:length(plot_z)
        field = strcat('Cooling_EnergyTransfer_Zone_',zonenames{plot_z(p)},'_J__Hourly_');
        field = field(1:min(63,length(field)));
        if any(strcmp(field,F))
            x = e_plus.(field)(index)/3600/1000;
            y = max(0,-Q_transfer(:,plot_z(p))/1000);
            y_lab = strcat(zonenames{plot_z(p)},' zone cooling (kW)');
            plot_fig(p_n+p,date2(2:end)-d1,x,y,y_lab)
        end
    end
    p_n = p_n + length(plot_z);
    
    for p = 1:1:length(plot_z)
        field = strcat('Heating_EnergyTransfer_Zone_',zonenames{plot_z(p)},'_J__Hourly_');
        field = field(1:min(63,length(field)));
        if any(strcmp(field,F))
            x = e_plus.(field)(index)/3600/1000;
            y = max(0,Q_transfer(:,plot_z(p))/1000);
            y_lab = strcat(zonenames{plot_z(p)},' zone heating (kW)');
            plot_fig(p_n+p,date2(2:end)-d1,x,y,y_lab)
        end
    end
    
    p_n = p_n + length(plot_z);
    for p = 1:1:length(plot_z)
        field = strcat(zonenames{plot_z(p)},'_ZoneAirTemperature_C__Hourly_');
        field = field(1:min(63,length(field)));
        if any(strcmp(field,F))
            x = e_plus.(field)(index);
            y = T_zone(2:end,plot_z(p));
            y_lab = strcat(zonenames{plot_z(p)},' zone temperature (C)');
            plot_fig(p_n+p,date2(2:end)-d1,x,y,y_lab)
        end
    end
    
    p_n = p_n + length(plot_z);
    for p = 1:1:length(plot_z)
        figure(p_n+p)
        field = strcat(zonenames{plot_z(p)},'_ZoneAirHeatBalanceSurfaceConvectionRate_W__Hourly_');
        field = field(1:min(63,length(field)));
        if any(strcmp(field,F))
            x = e_plus.(field)(index)/1000;
            y = Q_test_plots.surface_conv(:,plot_z(p))/1000;
            y_lab = strcat(zonenames{plot_z(p)},' zone surface convection gain (kW)');
            plot_fig(p_n+p,date2(2:end)-d1,x,y,y_lab)
        end
    end
   
%     %don't plot all internal gains so choose from list
%     g_list = {'LatentGain';'ConvectiveHeating';'RadiantHeating';'VisibleRadiationHeating';};
%     g_list2 = {'latent';'convected';'radiant';'visible';};
%     p_n = p_n + length(plot_z);
%     for k = 1:1:4
%         for p = 1:1:length(plot_z)
%             field = strcat(zonenames{plot_z(p)},'_ZoneTotalInternal',g_list{k},'Rate_W__Hourly_');
%             field = field(1:min(63,length(field)));
%             if any(strcmp(field,F))
%                 x = e_plus.(field)(index)/1000;
%                 y = loads.(g_list2{k})(:,plot_z(p))/1000;
%                 y_lab = strcat(zonenames{plot_z(p)},' zone ',g_list2{k},' gain (kW)');
%                 plot_fig(p_n+p,date2(2:end)-d1,x,y,y_lab)
%             end
%         end
%         p_n = p_n + length(plot_z);
%     end
%     
%     net_gain = zeros(length(index),length(zonenames));
%     for p = 1:1:length(plot_z)
%         field = strcat(zonenames{plot_z(p)},'_ZoneAirHeatBalanceOutdoorAirTransferRate_W__Hourly_');
%         field = field(1:min(63,length(field)));
%         if any(strcmp(field,F))
%             x = e_plus.(field)(index)/1000;
%             y = Q_test_plots.Infiltration(:,plot_z(p))/1000;
%             y_lab = strcat(zonenames{plot_z(p)},' zone outside air gain (kW)');
%             plot_fig(p_n+p,date2(2:end)-d1,x,y,y_lab)
%             net_gain = net_gain + x;
%         end
%     end
%     p_n = p_n + length(plot_z);
%     
%     for p = 1:1:length(plot_z)
%         field = strcat(zonenames{plot_z(p)},'_ZoneAirHeatBalanceInternalConvectiveHeatGainRate_W__Hour');
%         field = field(1:min(63,length(field)));
%         if any(strcmp(field,F))
%             x = e_plus.(field)(index)/1000;
%             y = gains.zone_sensible(:,plot_z(p))/1000;
%             y_lab = strcat(zonenames{plot_z(p)},' zone convection gain with ref, rack, water equip (kW)');
%             plot_fig(p_n+p,date2(2:end)-d1,x,y,y_lab)
%             net_gain(:,plot_z(p)) = net_gain(:,plot_z(p)) + x;
%         end
%     end
%     p_n = p_n + length(plot_z);
%     
%     for p = 1:1:length(plot_z)
%         y = Q_test_plots.net_internal_zone_gain/1000;
%         net_gain(:,plot_z(p)) = net_gain(:,plot_z(p)) + gains.ref_rack_water_latent(:,plot_z(p));
%         field = strcat(zonenames{plot_z(p)},'_ZoneTotalInternalLatentGainRate_W__Hourly_');
%         field = field(1:min(63,length(field)));
%         if any(strcmp(field,F))
%             x = e_plus.(field)(index)/1000;
%             net_gain(:,plot_z(p)) = net_gain(:,plot_z(p)) + x;
%             y_lab = strcat(zonenames{plot_z(p)},' net internal gain (kW)');
%             plot_fig(p_n+p,date2(2:end)-d1,net_gain(:,plot_z(p)),y,y_lab)
%         end
%     end
end
end%Ends function building_validate


function plot_fig(p,d,x,y,y_lab)
figure(p)
COD = 1 - sum((y - x).^2)/sum((x-mean(x)).^2);%Coefficient of determination
COD_B = 1 - sum((y - x).^2)/sum((y-mean(y)).^2);
plot(d,y,'k');
hold on
plot(d,x,'r');
legend({'WSU';'EnergyPlus';})
xlabel('Day of year')
ylabel(y_lab)
annotation('textbox',[.9 .75 .06 .1],'String',strcat('r^2 value :',num2str(max(COD,COD_B))),'EdgeColor','none')
annotation('textbox',[.9 .25 .06 .1],'String',strcat('NF value :',num2str(sum(y)/sum(x))),'EdgeColor','none')
end