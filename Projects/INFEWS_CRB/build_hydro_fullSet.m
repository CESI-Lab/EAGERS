function [Plant,TestData] = build_hydro_fullSet(name)
%% [Plant,TestData] = build_hydro_fullSet('FullColumbiaBasin');
%           All hydro Storage should be in kAcre-ft
%           All flows should be in Kcfs
%           All hydro power get converted from MW to kW
    
Plant.Name = name;    
%Setup Plant options
res = 1;%hourly resolution of test data
firstday = datenum([2007 1 1 1 0 0]); 
lastday = datenum([2008 10 1 0 0 0]);  
TestData.Timestamp= (firstday:res/24:lastday)';
TestData.Weather.Tdb = 20*ones(length(TestData.Timestamp),1);

%optimization options
Plant.optimoptions.Interval = 365; %Simulation Period
Plant.optimoptions.Horizon = 24*7; %Dispatch Horizon
Plant.optimoptions.Resolution = 1; %Initial Step Size
Plant.optimoptions.Topt = 600;
Plant.optimoptions.Tmpc = 60;
Plant.optimoptions.scaletime = 1;
Plant.optimoptions.tspacing = 'constant';
Plant.optimoptions.excessHeat = true;
Plant.optimoptions.excessCool = true;
Plant.optimoptions.thresholdSteps = 1;
Plant.optimoptions.method = 'Dispatch';
Plant.optimoptions.MixedInteger = true;
Plant.optimoptions.SpinReserve = false;
Plant.optimoptions.SpinReservePerc = 0;
Plant.optimoptions.forecast = 'Perfect';
Plant.optimoptions.solver = 'quadprog';
Plant.optimoptions.mode = 'virtual';
Plant.optimoptions.endSOC = 'Flexible';

%% Insert electric and gas utilities
dir=strrep(which('build_hydro_fullSet.m'),fullfile('Projects','INFEWS_CRB','build_hydro_fullSet.m'),'');
ut = 2;
load(fullfile(dir,'System_Library','Utility','Seattle1.mat'));
Plant.Generator(1,1).Type = 'Utility';
Plant.Generator(1,1).Name = 'Seattle_Light';
Plant.Generator(1,1).Source = 'Electricity';
Plant.Generator(1,1).Output = [];
Plant.Generator(1,1).Size = 0;
Plant.Generator(1,1).Enabled = 1;
Plant.Generator(1,1).VariableStruct = component;%.Name = 'Seattle_Light'; Plant.Generator(1,1).VariableStruct.SumStartMonth = 6; Plant.Generator(1,1).VariableStruct.SumStartDay = 1; Plant.Generator(1,1).VariableStruct.WinStartMonth = 10; Plant.Generator(1,1).VariableStruct.WinStartDay = 1;
Plant.Generator(1,1).VariableStruct.NodeName.Electric = 'Seattle';
Plant.Generator(1,1).VariableStruct.SellBackPerc = 50;% percent of purchase rate
Plant.Generator(1,1).VariableStruct.SellBackRate = 0.05;% $/kW or (-1) makes it sell back at purchase rate
Plant.Generator(1,1).VariableStruct.Latitude = 47.6062;
Plant.Generator(1,1).VariableStruct.Longitude =  -122.3321; 
Plant.Generator(1,1).VariableStruct.TimeZone =  -6;

Plant.Generator(2,1).Type = 'Utility';
Plant.Generator(2,1).Name = 'Northwest_Natural';
Plant.Generator(2,1).Source = 'NG';
Plant.Generator(2,1).Output = [];
Plant.Generator(2,1).Size = 0;
Plant.Generator(2,1).Enabled = 1;
Plant.Generator(2,1).VariableStruct.Timestamp = [datenum([2007 1 1 1 0 0]);datenum([2008 1 1 1 0 0])];%must be vector spanning 1 year (any year)
Plant.Generator(2,1).VariableStruct.Rate = [15;15];% $/MMBTU 

%% Plant Now at Seattle and Boise 
%% upgrade to read in electric generators from somewhere
elec_node_names = {'Seattle';'Boise';};
elec_station_names = {'Seattle Washington Electric Station';'Boise Idaho Electric Station';};
elec_station_long = [-122.3321;-116.2023;];
elec_station_lat = [47.6062;43.6150;];

for k = 1:1:length(elec_station_names)
    j = k+ut;
    Plant.Generator(j,1).Type = 'Electric Generator';
    Plant.Generator(j,1).Name = elec_station_names{k};
    Plant.Generator(j,1).Source = 'NG';
    Plant.Generator(j,1).Output.Capacity = (0:.1:1)';
    Plant.Generator(j,1).Output.Electricity = 0.6*ones(11,1);
    Plant.Generator(j,1).Size = 2e6;
    Plant.Generator(j,1).Enabled = 1;
    Plant.Generator(j,1).VariableStruct.NodeName.Electric = elec_node_names{k};
    Plant.Generator(j,1).VariableStruct.Longitude = elec_station_long(k); 
	Plant.Generator(j,1).VariableStruct.Latitude =  elec_station_lat(k);
    Plant.Generator(j,1).VariableStruct.TimeZone =  -6;
    Plant.Generator(j,1).VariableStruct.StateSpace.A = [0,1;-1.0823e-6,-.0021];Plant.Generator(j,1).VariableStruct.StateSpace.B = [0;1.0823e-6];Plant.Generator(j,1).VariableStruct.StateSpace.C = [1,0];Plant.Generator(j,1).VariableStruct.StateSpace.D = 0;
    Plant.Generator(j,1).VariableStruct.Startup.Time = [0,1000];Plant.Generator(j,1).VariableStruct.Startup.Electricity = [0,100];Plant.Generator(j,1).VariableStruct.Startup.Input = [0,166.67];
    Plant.Generator(j,1).VariableStruct.Shutdown.Time = [0,1000];Plant.Generator(j,1).VariableStruct.Shutdown.Electricity = [100,0];Plant.Generator(j,1).VariableStruct.Shutdown.Input = [166.67,0];
    Plant.Generator(j,1).VariableStruct.StartCost = 1e4;
    Plant.Generator(j,1).VariableStruct.dX_dt = 1e6;
end


%Read in information from Dams
instream = readtable('InstreamRequirements.xlsx');
Velocity = 1.5 * 0.000621371 * 3600; %conversion from meters/second to miles/second to miles/hour; 1 meter = 0.000621371 miles; 3600 seconds = 1 hour
fi = readtable('FullSetInformation.xlsx'); %Takes a while (>1second), but only loaded in Plant building
n_d = length(fi.Name); %Number of Dams in system
upriver = cell(n_d,1);
for i = 1:1:n_d
    j = k+ut+i;
    Plant.Generator(j,1).Type = 'Hydro Storage';
    Plant.Generator(j,1).Name = fi.Name{i};
    
    Plant.Generator(j,1).Source = 'Water';
    Plant.Generator(j,1).Output = [];
    Plant.Generator(j,1).Size = fi.MaxCapacity_kAcre_ft_(i);%This would be storage capacity in kilo-acre-feet
    Plant.Generator(j,1).Enabled = 1;
    
    Plant.Generator(j,1).VariableStruct.StartWYstate = 0.50; % ??? may need to be fixed?
    Plant.Generator(j,1).VariableStruct.Reservoir = fi.ReservoirName{i};
    Plant.Generator(j,1).VariableStruct.River = fi.River{i};
    Plant.Generator(j,1).VariableStruct.MaxGenCapacity = fi.Capacity_MW_(i)*1000; %power in kW
    Plant.Generator(j,1).VariableStruct.RampUp = fi.RampUp_MW_hr_(i)*1000; %power change in kW/hr
    Plant.Generator(j,1).VariableStruct.RampDown = fi.RampDown_MW_hr_(i)*1000; %power change in kW/hr
    Plant.Generator(j,1).VariableStruct.MaxGenFlow = fi.MaxGenFlow_kcfs_(i); %flow in kcfs
    Plant.Generator(j,1).VariableStruct.MaxSpillFlow = fi.MaxSpillFlow_kcfs_(i); %flow in kcfs 
    Plant.Generator(j,1).VariableStruct.MaxHead = fi.MaxHead_ft_(i); %feet
    Plant.Generator(j,1).VariableStruct.MinHead = fi.MinHead_ft_(i); %guess minimum height allowable
%     Plant.Generator(j,1).VariableStruct.Upriver = split_upstream(fi.Upstream{i});
    Plant.Generator(j,1).VariableStruct.Downriver = fi.Downstream{i};
    Plant.Generator(j,1).VariableStruct.Latitude = fi.LatDD(i);
    Plant.Generator(j,1).VariableStruct.Longitude = fi.LongDD(i);
    Plant.Generator(j,1).VariableStruct.TimeZone =  -6;
    Plant.Generator(j,1).VariableStruct.NodeName.Electric = fi.ElectricNode{i};
    Plant.Generator(j,1).VariableStruct.NodeName.Hydro = fi.VICName{i};
    if ~isempty(fi.InstreamConstraintSchedule{i})
        Plant.Generator(j,1).VariableStruct.InstreamConstraint = instream.(fi.InstreamConstraintSchedule{i});
    else 
        Plant.Generator(j,1).VariableStruct.InstreamConstraint = zeros(12,1);
    end
    Plant.Generator(j,1).VariableStruct.Time2Sea = fi.RiverMile_upstreamFromSea_(i)/Velocity; %miles/mph = hour
    
    if fi.IncludedInColSim_(i)==fi.inColSim_(i) %Check to see if included or confusion
        Plant.Generator(j,1).VariableStruct.ColSim.included = fi.inColSim_(i);
    else
        Plant.Generator(j,1).VariableStruct.ColSim.included = 'Error: Conflict';
    end
    Plant.Generator(j,1).VariableStruct.ColSim.naturalized_flow = fi.natFlow_(i);
    Plant.Generator(j,1).VariableStruct.ColSim.agriculture_withdrawl = fi.ColAgWithd(i);
    Plant.Generator(j,1).VariableStruct.ColSim.bias_corrected = fi.BiasCorrColSim_(i);
  
    Plant.Generator(j,1).VariableStruct.VIC.name = fi.VICName{i};
    Plant.Generator(j,1).VariableStruct.VIC.ID = fi.VICID(i);
    upriver(i) = {split_upstream(fi.Upstream{i})};
end


n_g = length(Plant.Generator);
e_nodes = {};
ee = 0;
for i = 1:1:n_g
    if isfield(Plant.Generator(i).VariableStruct,'NodeName') && isfield(Plant.Generator(i).VariableStruct.NodeName,'Electric')
        j = nonzeros((1:ee)'.*strcmpi(Plant.Generator(i).VariableStruct.NodeName.Electric,e_nodes));
        if isempty(j)
            ee = ee + 1;
            Plant.Network(ee).name = Plant.Generator(i).VariableStruct.NodeName.Electric;
            e_nodes(end+1,1) = {Plant.Network(ee).name};
            Plant.Network(ee).Equipment = {strcat(Plant.Generator(i).Type,'.',Plant.Generator(i).Name)};
            Plant.Network(ee).Electrical.connections = {};
            Plant.Network(ee).Electrical.Trans_Eff = [];
            Plant.Network(ee).Electrical.Trans_Limit = [];
            Plant.Network(ee).Electrical.Load = ee;
            Plant.Network(ee).Location.Longitude = Plant.Generator(i).VariableStruct.Longitude;
            Plant.Network(ee).Location.Latitude = Plant.Generator(i).VariableStruct.Latitude;
            Plant.Network(ee).Location.TimeZone = Plant.Generator(i).VariableStruct.TimeZone; 
        elseif ~strcmp(Plant.Generator(i).Type,'Hydro Storage')
            Plant.Network(j).Equipment = [Plant.Network(j).Equipment;strcat(Plant.Generator(i).Type,'.',Plant.Generator(i).Name)];
        end
    end
end
Plant.Network(1).Electrical.connections = {'Boise'};
Plant.Network(1).Electrical.Trans_Eff = 0.999;
Plant.Network(1).Electrical.Trans_Limit = inf;
            
h_nodes = {};
hh = ee;
for i = 1:1:n_g
    if isfield(Plant.Generator(i).VariableStruct,'NodeName') && isfield(Plant.Generator(i).VariableStruct.NodeName,'Hydro')
        j = nonzeros((ee+1:hh)'.*strcmpi(Plant.Generator(i).VariableStruct.NodeName.Hydro,h_nodes));
        if isempty(j)
            hh = hh + 1;
            Plant.Network(hh).name = Plant.Generator(i).VariableStruct.NodeName.Hydro;
            h_nodes(end+1,1) = {Plant.Network(hh).name};
            Plant.Network(hh).Equipment = {strcat(Plant.Generator(i).Type,'.',Plant.Generator(i).Name)};
            Plant.Network(hh).Electrical.connections = {Plant.Generator(i).VariableStruct.NodeName.Electric};
            Plant.Network(hh).Electrical.Trans_Eff = 1;
            Plant.Network(hh).Electrical.Trans_Limit = inf;
            Plant.Network(hh).Electrical.Load = [];
            Plant.Network(hh).Hydro.connections = {};
            Plant.Network(hh).Hydro.Time2Sea = Plant.Generator(i).VariableStruct.Time2Sea;
            Plant.Network(hh).Hydro.InstreamFlow = Plant.Generator(i).VariableStruct.InstreamConstraint;
            Plant.Network(hh).Location.Longitude = Plant.Generator(i).VariableStruct.Longitude;
            Plant.Network(hh).Location.Latitude = Plant.Generator(i).VariableStruct.Latitude;
            Plant.Network(hh).Location.TimeZone = Plant.Generator(i).VariableStruct.TimeZone; 
        else
            Plant.Network(j).Equipment = [Plant.Network(j).Equipment;strcat(Plant.Generator(i).Type,'.',Plant.Generator(i).Name)];
        end
    end
end
%% need second loop so that downstream node is created, so as to find it
for i = 1:1:n_g
    if isfield(Plant.Generator(i).VariableStruct,'NodeName') && isfield(Plant.Generator(i).VariableStruct.NodeName,'Hydro')
        j = nonzeros((ee+1:hh)'.*strcmpi(Plant.Generator(i).VariableStruct.NodeName.Hydro,h_nodes));
        if ~isempty(Plant.Generator(i).VariableStruct.Downriver)%find downstream node
            k = nonzeros((ee+1:hh)'.*strcmpi(Plant.Generator(i).VariableStruct.Downriver,h_nodes));%check node names
            if isempty(k)%check equipment names
                for k = ee+1:hh
                    if any(strcmp(strcat('Hydro Storage.',Plant.Generator(i).VariableStruct.Downriver),Plant.Network(k).Equipment))
                        break
                    end
                end
            end
            Plant.Network(j).Hydro.connections = h_nodes(k-ee);%downstream node
        end
    end
end
%% make sure electrical connections are both ways
nodes = length(Plant.Network);
n_names = cell(nodes,1);
for i = 1:1:nodes
    n_names(i) = {Plant.Network(i).name};
end
for i = 1:1:nodes
    connect = Plant.Network(i).Electrical.connections;
    for j = 1:1:length(connect)
        k = nonzeros((1:nodes)'.*strcmpi(connect{j},n_names));
        if ~any(strcmp(n_names{i},Plant.Network(k).Electrical.connections))
            Plant.Network(k).Electrical.connections = [Plant.Network(k).Electrical.connections;n_names{i}];
            Plant.Network(k).Electrical.Trans_Eff = [Plant.Network(k).Electrical.Trans_Eff,Plant.Network(i).Electrical.Trans_Eff(j)];
            Plant.Network(k).Electrical.Trans_Limit = [Plant.Network(k).Electrical.Trans_Limit,Plant.Network(i).Electrical.Trans_Limit(j)];
        end
    end
end

%% Data input file
acre_ft_per_week_2_kcfs = 7.20238e-5;
nat = true;
if nat
    fd = readtable('vic_data.xlsx','Sheet','naturalized','Range','B2:AY29315');
    vic_date = readtable('vic_data.xlsx','Sheet','naturalized','Range','A2:A29315');
    [~,fd_river] = xlsread('vic_data.xlsx','naturalized','B1:AY1');
else
    fd = readtable('vic_data.xlsx','Sheet','bias_corrected','Range','B2:AY1615');
    vic_date = readtable('vic_data.xlsx','Sheet','bias_corrected','Range','A2:A1615');
    [~,fd_river] = xlsread('vic_data.xlsx','bias_corrected','B1:AY1');
end
vic_date = datenum(vic_date.Date_);
fd_names = cell(n_d,1);
old_names = fd.Properties.VariableDescriptions';
table_names = fd.Properties.VariableNames;
for i = 1:1:n_d
    a = strfind(old_names{i},'''');
    fd_names(i) = {old_names{i}(a(1)+1:a(2)-1)};
end

%% eliminate nan values
for i = 1:1:n_d
    if isnan(fd.(table_names{i})(1))
        k = nonzeros((1:n_d)'.*strcmp(fd_names{i},fi.Name));
        if ~isempty(upriver{k})
            for j = 1:1:length(upriver{k})
                up = nonzeros((1:n_d)'.*strcmp(upriver{k}{j},fd_names));
                if ~isnan(fd.(table_names{up})(1))
                    fd.(table_names{i}) = fd.(table_names{up});%replace nan with downstream flow
                end
            end
            if isnan(fd.(table_names{i})(1))%look further upstream
                up = nonzeros((1:n_d)'.*strcmp(upriver{up}{1},fd_names));
                if ~isnan(fd.(table_names{up})(1))
                    fd.(table_names{i}) = fd.(table_names{up});%replace nan with downstream flow
                else
                    down = nonzeros((1:n_d)'.*strcmp(fi.Downstream{k},fd_names));
                    if isempty(fi.Downstream{k}) || isnan(fd.(table_names{down})(1))%no downstream node, or downstream is nan as well
                        %take average of any other measurements on that river
                        same_r = nonzeros((1:n_d)'.*strcmp(fd_river(i),fd_river'));
                        disp('both upriver and downriver nodes are nan')
                    else
                        fd.(table_names{i}) = fd.(table_names{down});%replace nan with downstream flow
                    end
                end
            end
        else
            down = nonzeros((1:n_d)'.*strcmp(fi.Downstream{k},fd_names));
            fd.(table_names{i}) = fd.(table_names{down});%replace nan with downstream flow
        end
    end
end
%% scale flow
for i = 1:1:n_d
    fd.(table_names{i}) = fd.(table_names{i})*acre_ft_per_week_2_kcfs;
end

%Need SourceSink, InFlow
SS = zeros(length(vic_date),n_d);
InF = zeros(length(vic_date),n_d);
TestData.Hydro.SourceSink = zeros(length(TestData.Timestamp),n_d);
TestData.Hydro.InFlow = zeros(length(TestData.Timestamp),n_d);
for i = 1:1:n_d
    k = nonzeros((1:n_d)'.*strcmp(fi.Name{i},fd_names));
    InF(:,i) = fd.(table_names{k});
    SS(:,i) = fd.(table_names{k});
    for j = 1:1:length(upriver{i})
        up = nonzeros((1:n_d)'.*strcmp(upriver{i}{j},fd_names));
        SS(:,i) = SS(:,i) - fd.(table_names{up});%remove flow at upstream node
    end
    if mean(SS(:,i))<-.05*mean(fd.(table_names{k}))
        disp('potential error in source/sink calculation')
    end
%     TestData.Hydro.SourceSink(:,i) = interp1(vic_date,SS(:,i),TestData.Timestamp);
%     TestData.Hydro.InFlow(:,i) = interp1(vic_date,InF(:,i),TestData.Timestamp);
end
%%if test data frequency is higher than water data uses prior point, if lower frequency, average the points
TestData.Hydro.SourceSink = interp_data(vic_date,SS,TestData.Timestamp);
TestData.Hydro.InFlow = interp_data(vic_date,InF,TestData.Timestamp);

%% Demand
d_s = datevec(TestData.Timestamp(1));
year_1 = d_s(1);
d_s = datevec(TestData.Timestamp(end));
year_end = d_s(1);
years = year_1:1:year_end;
bpa_date = [];
total_load = [];
for i = 1:1:length(years)
    if years(i)<=2010
        bpa1 = readtable(strcat(num2str(years(i)),'_BPA_TotalLoad_5Min.xlsx'),'Range','A17:H55000');
        bpa2 = readtable(strcat(num2str(years(i)),'_BPA_TotalLoad_5Min.xlsx'),'Range','J17:Q55000');
        bpa_date = [bpa_date;datenum([bpa1.DateTime(1:find(isnan(bpa1.Month),1,'first')-1);bpa2.DateTime(1:find(isnan(bpa2.Month),1,'first')-1)])];
        total_load = [total_load;bpa1.TOTALBPACONTROLAREALOAD_MW_SCADA45583_(1:find(isnan(bpa1.Month),1,'first')-1);bpa2.TOTALBPACONTROLAREALOAD_MW_SCADA45583_(1:find(isnan(bpa2.Month),1,'first')-1)]; 
    else
        bpa1 = readtable(strcat(num2str(years(i)),'_BPA_TotalLoad_5Min.xlsx'),'Sheet','January-June','Range','A19:G55000');
        bpa2 = readtable(strcat(num2str(years(i)),'_BPA_TotalLoad_5Min.xlsx'),'Sheet','July-December','Range','A19:G55000');
        bpa_date = [bpa_date;datenum([bpa1.DateTime(1:find(isnan(bpa1.Month),1,'first')-1);bpa2.DateTime(1:find(isnan(bpa2.Month),1,'first')-1)])];
        total_load = [total_load;bpa1.TOTALBPACONTROLAREALOAD_MW_SCADA45583_(1:find(isnan(bpa1.Month),1,'first')-1);bpa2.TOTALBPACONTROLAREALOAD_MW_SCADA45583_(1:find(isnan(bpa2.Month),1,'first')-1)]; 
    end
end
total_load = total_load*1000; %convert MW to kW
TestData.Demand.E(:,1) = interp_data(bpa_date,0.5*total_load,TestData.Timestamp);
TestData.Demand.E(:,2) = TestData.Demand.E(:,1);
Plant.Network(1).Electrical.Load = 1; %Seattle
Plant.Network(2).Electrical.Load = 2; %Boise
end

function select_data = interp_data(timestamp,data,date)
%% code taken from get_data, but modified to take prior point instead of interpolate when date timestep is smaller than data timestep
n_s = length(date);
index = nnz(timestamp<date(1)+1e-6)+ (0:n_s-1)'; %+1e-6 to avoid rounding problems
need_interp = true;
if length(timestamp)>=index(end)
    r = (timestamp(index(2:end))-timestamp(index((1:end-1))))./(date(2:end)-date(1:end-1));
    if all(abs(r-1)<1e-7)
        need_interp = false;
    end
end
 if ~need_interp%no need to interpolate
    select_data = data(index,:);
else
    select_data = zeros(n_s,length(data(1,:)));
    x1 = index(1);
    for t = 1:1:n_s
        k = 0;
        r0 = (timestamp(x1+1) - date(t))/(timestamp(x1+1) - timestamp(x1));
        while r0<0 %summing multiple data points
            k = k + 1+r0;
            select_data(t,:) = select_data(t,:) + (1+r0)*data(x1,:);
            x1 = x1+1;
            r0 = (timestamp(x1+1) - date(t))/(timestamp(x1+1) - timestamp(x1));
        end
        if k>0%summing multiple points                
            select_data(t,:) = (select_data(t,:) + (1-r0)*data(x1,:))/(k+(1-r0));
        else %using prior point
            select_data(t,:) =  data(x1,:);
%         else %interpolating
%             select_data(t,:) =  (r0*data(x1,:) + (1-r0)*data(x1+1,:));
        end
    end
end
end%Ends function interp_data

function upriver = split_upstream(text)
if strcmp(text,'Source')
    upriver = {};
else
    com = strfind(text,',');
    if isempty(com)
        upriver = {text};
    else
        upriver = cell(length(com)+1,1);
        upriver(1) = {text(1:com(1)-1)};
        for j = 2:length(com)
            a = text(com(end)+1:end); 
            if strcmp(a(1),' ')
                a = a(2:end);
            end
            upriver(j) = {a};
        end
        a = text(com(end)+1:end); 
        if strcmp(a(1),' ')
            a = a(2:end);
        end
        upriver(end) = {a};
    end
end
end%Ends function