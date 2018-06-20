%% comparison for paper
global Plant TestData
load('CampusRan_2_10_18.mat');
MIPlant = Plant;
TestData = [];
load_test_data
TestData.RealTimeData = interpolate_data(TestData,Plant.optimoptions.Resolution*3600,0.00);%create test data at correct frequency
%% Run non-mixed integer version from same IC
Plant.optimoptions.MixedInteger = false;
NumSteps = nnz(Plant.Dispatch.Timestamp)-1;
timers = zeros(NumSteps,3); % To record times set to zeros(1,3), to not record set to [];
nG = length(Plant.Generator);
LB = zeros(1,nG);
for i=1:1:nG
    if Plant.OpMatA.Organize.Dispatchable(i)
        states = Plant.Generator(i).QPform.states(1:nnz(~cellfun('isempty',Plant.Generator(i).QPform.states(:,end))),end);
        for j = 1:1:length(states)
            LB(i) = LB(i) + Plant.Generator(i).QPform.(states{j}).lb(2);
        end
    end
end
if isfield(Plant,'Building') && ~isempty(Plant.Building)
    Buildings = Plant.Building;
else
    Buildings = [];
end 
Time = build_time_vector(Plant.optimoptions);%% set up vector of time interval
Solution.Dispatch = [];
Si = 1;
Date = Plant.Dispatch.Timestamp(1)+[0;Time/24];
DispatchWaitbar=waitbar(0,'Running Dispatch','Visible','on');
while Si<NumSteps-1
    for i = 1:1:nG
        Plant.Generator(i).CurrentState = MIPlant.Dispatch.GeneratorState(Si,i);   
        Plant.Generator(i).Status = Plant.Generator(i).CurrentState>LB(i);
    end
    [Forecast,Plant.Generator,Buildings] = update_forecast(Plant.Generator,Buildings,Plant.subNet,Plant.optimoptions,Date(2:end));%% function that creates demand vector with time intervals coresponding to those selected
    [Solution,~] = dispatch_loop(Plant.Generator,Buildings,Plant.subNet,Plant.OpMatA,Plant.OpMatB,Plant.OneStep,Plant.optimoptions,Date,Forecast,Solution);
    timers(Si,:) = Solution.timers;
    [C,~,~] = net_cost(Plant.Generator,Solution.Dispatch,Date,'Dispatch');
    Plant.Predicted.Cost(Si) = sum(C);
    [Si,Date,Plant.Generator,Plant.Building,Hydro] = dispatch_record(Plant.Generator,Buildings,Plant.subNet,Plant.optimoptions,TestData.RealTimeData,Si,Date,Forecast,Solution);
    if ~isempty(Hydro)
        nS = length(Hydro(:,1));
        TestData.RealTimeData.Hydro.OutFlow(Si-nS+1:Si,:) = Hydro;
    end
    waitbar(Si/NumSteps,DispatchWaitbar,strcat('Running Dispatch'));
end
cQPlant = Plant;

MIDisp = MIPlant.Plant.Dispatch.GeneratorState(2:NumSteps,:);
cQPDisp = cQPlant.Plant.Dispatch.GeneratorState(2:NumSteps,:);
ICMI = MIPlant.Plant.Dispatch.GeneratorState(1,:);
ICcQP = cQPlant.Plant.Dispatch.GeneratorState(1,:);

%find difference in dispatch setpoints
dispAble = 1:16;
ndisp = length(dispAble);
samesame = (MIDisp<cQPDisp+1) & (MIDisp>cQPDisp-1);
percSame = nnz(samesame(:,dispAble))/(NumSteps*ndisp);

%find differences in dispatch on/off
gensOnOff = 3:13;%units that have non-zero lb
ngens = length(gensOnOff);
MIonoff = (MIDisp(:,gensOnOff)>0);
cQPonoff = (cQPDisp(:,gensOnOff)>0);

sameonoff = (MIonoff==cQPonoff);
percOnOff = nnz(sameonoff)/(NumSteps*ngens);

%plot dispatches
% for i = 1:1:length(dispAble)
%     figure(i)
%     plot(MIDisp(:,dispAble(i)));
%     hold on
%     plot(cQPDisp(:,dispAble(i)));
%     title(strcat('disp for gen ',PlantcQP.Plant.Generator(dispAble(i)).Name));
%     legend('MI','cQP')
%     hold off
% end

%find difference in cost
DemandE = MIPlant.Plant.Dispatch.Demand.E(2:NumSteps);
Timestamp = MIPlant.Plant.Dispatch.Timestamp(2:NumSteps);
method = 'Dispatch';
[MICost,MImissedDem] = net_cost(Plant.Generator,MIDisp,Timestamp,method);
[cQPCost,cQPmissedDem] = net_cost(Plant.Generator,cQPDisp,Timestamp,method);
%cost difference
costDiff = cQPCost-MICost;
muCost = mean(costDiff);
sigmaCost = std(costDiff);
%individual method costs
mucQPCost = mean(cQPCost);
sigmacQPCost = std(cQPCost);
muMICost = mean(MICost);
sigmaMICost = std(MICost);
%missed demand due to difference between linear and non-linear chiller
%efficiency
muMIDem = mean(MImissedDem);
sigmaMIDem = std(MImissedDem);
mucQPDem = mean(cQPmissedDem);
sigmacQPDem = std(cQPmissedDem);

%plot 
%create gaussian of costs
nsigmas = 4;
dsigma = 1/8;
nbins = nsigmas*2*1/dsigma+1;

cQPcostBins = zeros(nbins*2-2,1);%bin by sigma/2, go out to six sigma +1 in each direction
MIcostBins = zeros(nbins*2-2,1);
%find bottom of bins
lbBincQP = mucQPCost-sigmacQPCost*dsigma/2-sigmacQPCost*nsigmas+[0:1:nbins-1]*sigmacQPCost*dsigma;
lbBinMI = muMICost-sigmaMICost*dsigma/2-sigmaMICost*nsigmas+[0:1:nbins-1]*sigmaMICost*dsigma;
%find top of bins
ubBincQP = mucQPCost+sigmacQPCost*dsigma/2-sigmacQPCost*nsigmas+[0:1:nbins-1]*sigmacQPCost*dsigma;
ubBinMI = muMICost+sigmaMICost*dsigma/2-sigmaMICost*nsigmas+[0:1:nbins-1]*sigmaMICost*dsigma;
for i = 1:1:nbins%count number in bins
    %find number in bins
    cQPcostBins((i-1)*2+1:i*2) = nnz(and((cQPCost>lbBincQP(i)),(cQPCost<=ubBincQP(i))));
    MIcostBins((i-1)*2+1:i*2) = nnz(and((MICost>lbBinMI(i)),(MICost<=ubBinMI(i))));
    %lbBincQP = lbBincQP+sigmacQPCost/2;%move up one sigma
    %ubBincQP = ubBincQP+sigmacQPCost/2;
end
%plot gaussian of cost
figure(1)
costscQP = sort([lbBincQP,ubBincQP]);
costsMI = sort([lbBinMI,ubBinMI]);
plot(costsMI, MIcostBins);
hold on
plot(costscQP,cQPcostBins);
title('Cost of Distribution of Dispatches')
xlabel('Cost in $ at each hour')
ylabel('number of hours in a year')
legend('mcQP', 'cQP');
xlim([0,2500]);
hold off    

%plot cost of one day
iday = 4345;
figure(2)
plot(MICost(iday:iday+24))
hold on
plot(cQPCost(iday:iday+24))
hold off
title('Hourly Cost of Actual Dispatch for July 1st')
xlabel('Time in hours')
ylabel('Cost of disptach in $')
legend('mcQP', 'cQP')
JulyMICost = sum(MICost(iday:iday+24))
JulycQPCost = sum(cQPCost(iday:iday+24))
%show that different costs are related to different on/off, but that total
%cost is similar
Time = [1:24];
%plotComparison(MIDisp(iday:iday+24,:),Plant.Dispatch.Demand,Time,3)
%plotComparison(cQPDisp(iday:iday+24,:),Plant.Dispatch.Demand,Time,6)


% %find out when fuel cells are at less than 99.5% of max power
% FC1dwnMI = nnz(MIDisp(:,5)<1990);
% FC2dwnMI = nnz(MIDisp(:,6)<1990);
% FC1dwncQP = nnz(cQPDisp(:,5)<1990);
% FC2dwncQP = nnz(cQPDisp(:,6)<1990);
% FConoff = nnz(MIonoff(:,[3,4])==cQPonoff(:,[3,4]));
% FCsame = nnz((MIDisp(:,[5,6])<cQPDisp(:,[5,6])+1).*(MIDisp(:,[5,6])>cQPDisp(:,[5,6])-1));
% %percents
% FConPerc = FConoff/(endtime*2);
% FCsamePerc = FCsame/(endtime*2);

%find distance between generators of like types
%gas turbines 


%% compare concurrent run
%since the non-concurrent version shows different peaks, need to analyze
%horizon dispatches
%plot predicted dispatch costs
horizonCMI = MIPlant.Plant.Predicted.Cost(1:NumSteps-1);
horizonCcQP = cQPlant.Plant.Predicted.Cost(1:NumSteps-1);
horizonCcQP_noIC = zeros(NumSteps-1,1);
for i = 1:1:NumSteps-1
    horizonCcQP_noIC(i) = sum(net_cost(Plant.Generator,cQPlant.Plant.Predicted.GenDisp(:,:,i),cQPlant.Plant.Predicted.Timestamp(:,i),method));
end
figure(9)
plot(horizonCMI)
hold on
plot(horizonCcQP)
hold off
title('Cost of Horizon Dispatch')
ylabel('Cost of Dispatch for 24 hour Horizon in $')
xlabel('Hour of the year')
xlim([0,8714])
legend('mcQP', 'cQP')
muhorizonCMI = mean(horizonCMI);
sigmahCMI = std(horizonCMI);
muhorizonCcQP = mean(horizonCcQP);
sigmahCcQP = std(horizonCcQP);

muCcQP_CMI = mean(horizonCMI-horizonCcQP);
muCcQPnoIC_CMI = mean(horizonCMI-horizonCcQP_noIC);

%compare starts
MIstarts = MIonoff(2:end,:)>MIonoff(1:end-1,:);
cQPstarts = cQPonoff(2:end,:)>cQPonoff(1:end-1,:);
nMIstarts = nnz(MIstarts);
ncQPstarts = nnz(cQPstarts);




