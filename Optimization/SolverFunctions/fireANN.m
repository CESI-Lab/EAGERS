function Locked = fireANN(PrevDisp,Forecast,scaleCost,training)
%%
% fireANN is an upgrade of the previous fireNeurons. fireNeurons only
% conducted unit commitment, so one of the features was the dispatch from
% fitA. fireANN goes directly from the demand update step to the unit
% commitment, so features are the previous timesteps' dispatch, the
% forecasted demand, the initial conditions, and the scale costs. This
% function creates and trains an ANN for unit commitment and then uses that
% trained ANN for dispatching in later timesteps. Contact author for
% fireNuerons (in use before 2/20/2018)
%
% fireANN(previousDispatch, Forecasted demand, scalecosts, training)
% 
% INPUTS:
%   PrevDisp    -   The dispatch from the previous timestep in matrix form
%                   with size (#timesteps X #components)
%   Forecast    -   The forecasted demand in a structure where Forecast.Demand.E
%                   is the electrical demand in (#timesteps X 1),
%                   Forecast.Demand.H is the heat demand in (#timesteps X 1), and
%                   Forecast.Demand.C is the cooling demand in (#timesteps X 1)
%   scaleCost   -   This is the scaling for the cost of power for each
%                   component for each timestep. It is proportional to the
%                   cost of the utility and the cost of fuel
%   training    -   This is a boolean value denoting if you need to
%                   initialize and train a neural network, or if you have   
%                   already created and trained one.
%
% OUTPUTS:
%   Locked      -   Matrix of ones and zeros denoting online or offline for
%                   each component at each timestep size (#timesteps X
%                   #components)
%
% Last Updated 2/21/2018 by Nadia Panossian
%%

global Plant Model_dir TestData
global lockedNet lockederror %these are created in this function and used in each subsequent step

%% initialize crucial values
Time = build_time_vector(Plant.optimoptions);%% set up dt vector of time interval length
dt = [Time - [0; Time(1:end-1)]]';
nG = length(Plant.Generator);
nS = length(Time);
GenDemand  = Forecast.Demand;%only need the demand portion
Out = fieldnames(GenDemand); %type of energy outputs
renews = []; %number of renewables
gens = []; %generators
stor = []; %storage
util = [];
UB = zeros(1,nG);
LB = zeros(1,nG);
QPsingle = Plant.OpMatB;
for i = 1:1:nG%find renewables and upper and lower bounds
    states = QPsingle.Organize.States{i};
    if strcmp(Plant.Generator(i).Source, 'Renewable')
        renews(end+1) = i;
        UB(i) = Plant.Generator(i).Size;
        LB(i) = 0;
    elseif ~isempty(strfind(Plant.Generator(i).Type,'Storage'))
        stor(end+1) = i;
        UB(i) = Plant.Generator(i).Size;
        LB(i) = UB(i).*.2;
    elseif ~strcmp(Plant.Generator(i).Type,'Utility') %&& ~strcmp(Plant.Generator(i).Type,'Heater')%no unit commitment for utilities, storage, renewables, or heaters
        gens(end+1) = i;
        UB(i) = Plant.Generator(i).Size;
        LB(i) = min(0,sum(QPsingle.lb(states)));
    else
        UB(i) = inf;
        LB(i) = 0;
        util(end+1) = i;
    end
end 
ngens = length(gens);
%specify number of neural network input features
nfeatures = 3*ngens+1+length(stor)+length(Out);

%% create training data for 365 days
if training %train the first time through
    %select which file to train from
    filename = SelectTrainingData;
    %select 2_10_18
    if ~isempty(filename)
        PlantTrain = load(fullfile(Model_dir,filename));        
    else
        %save current conditions
        Gen = Plant.Generator;
        Build = Plant.Building;
        Plant.optimoptions.solver = 'quadprog';
        
        %run through a year for training
        disp('Running Simulation for Training Examples')
        RunSimulation(TestData.Timestamp(1),[])
        PlantTrain = Plant;
        
        %reset values back to pre-training vals
        Plant.optimoptions.solver = 'ANN';
        Plant.optimoptions.forecast = forecastMethod;
        Plant.Generator = Gen;
        Plant.Building = Build;
    end
    
    if isfield(PlantTrain, 'Plant')
        PlantTrain = PlantTrain.Plant;
        %sort solutions into features
        disp('sorting training data')
        dt = dt';
        %initialize space for training info
        ndisps = length(PlantTrain.Predicted.GenDisp(1,1,:))-2;
        GenDisp = zeros(ndisps*nS,nG);
        scaleCost = zeros(ndisps*nS,1);
        RenewE = zeros(ndisps*nS,1);
        PrevDisp = zeros(ndisps*nS,nG);
        for j = 1:1:length(Out)
            GenDemand.(Out{j}) = zeros(ndisps*nS,1);
        end

        %load generator dispatch
        for t = 1:1:ndisps
            %don't include initial conditions, because they are covered
            %in the previous dispatch
            GenDisp((t-1)*nS+1:t*nS,:) = PlantTrain.Predicted.GenDisp(:,:,t);
        end

        %load utility cost
        sumTable = PlantTrain.Generator(util(1)).VariableStruct.SumRateTable;
        winTable = PlantTrain.Generator(util(1)).VariableStruct.WinRateTable;
        winRates = PlantTrain.Generator(util(1)).VariableStruct.WinRates;
        sumRates = PlantTrain.Generator(util(1)).VariableStruct.SumRates;
        Date0 = datevec(PlantTrain.Dispatch.Timestamp(1));
        year = Date0(1);
        sprgDate = datenum([year, PlantTrain.Generator(util(1)).VariableStruct.SumStartMonth, PlantTrain.Generator(util(1)).VariableStruct.SumStartDay]);
        fallDate = datenum([year, PlantTrain.Generator(util(1)).VariableStruct.WinStartMonth, PlantTrain.Generator(util(1)).VariableStruct.WinStartDay]);
        for t = 1:1:ndisps
            Timestamp = PlantTrain.Dispatch.Timestamp(t);
            day = weekday(Timestamp);
            if Timestamp<sprgDate || Timestamp>=fallDate
                scaleCost((t-1)*nS+1:t*nS) = (winRates(winTable(day,:))'.*dt);
            else
                scaleCost((t-1)*nS+1:t*nS) = (sumRates(sumTable(day,:))'.*dt);
            end
            %load renewable energy and electrical demand
            RenewE((t-1)*nS+1:t*nS) = sum(PlantTrain.Predicted.GenDisp(:,renews,t),2);
            for j = 1:1:length(Out)
                GenDemand.(Out{j})((t-1)*nS+1:t*nS) = PlantTrain.Data.Demand.(Out{j})((t-1)*4+1:4:(t+nS-1)*4);%convert from 15 min to hourly
            end
            %load the previous dispatch which is what was predicted a timestep ago
            if t>1
               PrevDisp((t-1)*nS+1:t*nS,:) = PlantTrain.Predicted.GenDisp(:,:,t-1);%the predictions don't include the gas utility
            else
               PrevDisp(1:nS,:) = PlantTrain.Predicted.GenDisp(:,:,end);%runs should start where they end (i.e. from Jan 1st to Dec 31st or from Monday to Sunday)
            end
        end
        %convert generator dispatches to OnOff
        OnOffT = (GenDisp(2:end,gens)>0);
    
    else%if you are using a pretrained ANN
        lockedNet = PlantTrain.lockedNet;
        training = false; %don't retrain
        ndisps = 1;
        scaleCost = scaleCost(:,util(1));
        RenewE = zeros(nS,1);%first step doesn't have renewables
        PrevDisp = ones(nS+1,nG);
    end
%if you are not training    
else
    ndisps = 1;
    scaleCost = scaleCost(:,util(1));
    RenewE = sum(PrevDisp(:,renews),2);
end

%% create, train, and dispatch a locked network
%put all features into input matrices
lockInputs = zeros(nS*ndisps, nfeatures);
parfor i = 1:1:length(Out)
    if strcmp(Out,'E')
        lockInputs(:,i) = (GenDemand.(Out{i})-RenewE)./(max(GenDemand.(Out{i})));
    else
        lockInputs(:,i) = GenDemand.(Out{i})./(max(GenDemand.(Out{i})));
    end
end
col = length(Out)+1;%collumn for input
%add cost terms
lockInputs(:,col) = scaleCost./max(max(Plant.Generator(util(1)).VariableStruct.SumRates));%cost of utility electrical
%add power from storage given fit A prediction
maxSPower = (ones(nS,1) * max(PrevDisp(1:end-1,stor)));%max(PrevDisp(1:end-1,stor)-PrevDisp(2:end,stor)));
for t = 1:1:ndisps
    %add power terms from storage prediction
    if training
        if t>1
            lockInputs((t-1)*nS+1, col+1:col+length(stor)) = (PrevDisp((t-2)*nS+1,stor)-PrevDisp((t-1)*nS+1,stor)) ./ maxSPower(1,:);
        else %use the initial condition for the first timestep, assume circular
            lockInputs(1,col+1:col+length(stor)) = (GenDisp(1,stor)-PrevDisp(1,stor)) ./ maxSPower(1,:);
        end
        lockInputs((t-1)*nS+2:t*nS,col+1:col+length(stor)) = (PrevDisp((t-1)*nS+1:t*nS-1,stor)-PrevDisp((t-1)*nS+2:t*nS,stor)) ./ maxSPower(2:end,:);%SOC of battery at previous step
    else
        lockInputs(1:nS,col+1:col+length(stor)) = (PrevDisp(1:end-1,stor)-PrevDisp(2:end,stor)) ./maxSPower;
    end
end

col = col+length(stor)+1;
for i = 1:1:ngens
    j = gens(i);
    j1 = QPsingle.Organize.States{j};
    if ~isempty(j1)
        %ramp rate terms are not important so are not included
        for t = 1:1:ndisps
                lockInputs((t-1)*nS+1:t*nS,col) = [PrevDisp((t-1)*nS+2:t*nS,j);PrevDisp((t-1)*nS+1,j)]./UB(j);%Prev dispatch estimate of now
                lockInputs((t-1)*nS+1:t*nS,col+1) = PrevDisp((t-1)*nS+1:t*nS,j)./UB(j);%Initial condition
                lockInputs((t-1)*nS+1:t*nS,col+2) = [PrevDisp((t-1)*nS+3:t*nS,j);PrevDisp((t-1)*nS+(1:2),j)]./UB(j);%Previous dispatch estimate of t+1
%                 lockInputs((t-1)*nS+1:t*nS,col+3) = [PrevDisp((t-1)*(nS+1)+4:t*(nS+1),j);PrevDisp((t-1)*(nS+1)+1:(t-1)*(nS+1)+2,j)]./UB(j);%Prev dispatch at t+2
%                 lockInputs((t-1)*nS+1:t*nS,col+4) = [PrevDisp((t-1)*(nS+1)+5:t*(nS+1),j);PrevDisp((t-1)*(nS+1)+1:(t-1)*(nS+1)+3,j)]./UB(j);%Prev dispatch at t+3
%                 lockInputs((t-1)*nS+1:t*nS,col+5) = [PrevDisp((t-1)*(nS+1)+6:t*(nS+1),j);PrevDisp((t-1)*(nS+1)+1:(t-1)*(nS+1)+4,j)]./UB(j);%Prev dispatch at t+4
%                 lockInputs((t-1)*nS+1:t*nS,col+6) = [PrevDisp((t-1)*(nS+1)+7:t*(nS+1),j);PrevDisp((t-1)*(nS+1)+1:(t-1)*(nS+1)+5,j)]./UB(j);%Prev dispatch at t+5
%                 lockInputs((t-1)*nS+1:t*nS,col+7) = [PrevDisp((t-1)*(nS+1)+8:t*(nS+1),j);PrevDisp((t-1)*(nS+1)+1:(t-1)*(nS+1)+6,j)]./UB(j);%Prev dispatch at t+6
%                 lockInputs((t-1)*nS+1:t*nS,col+8) = [PrevDisp((t-1)*(nS+1)+9:t*(nS+1),j);PrevDisp((t-1)*(nS+1)+1:(t-1)*(nS+1)+7,j)]./UB(j);%Prev dispatch at t+7
%                 lockInputs((t-1)*nS+1:t*nS,col+9) = [PrevDisp((t-1)*(nS+1)+10:t*(nS+1),j);PrevDisp((t-1)*(nS+1)+1:(t-1)*(nS+1)+8,j)]./UB(j);%Prev dispatch at t+8
%                 lockInputs((t-1)*nS+1:t*nS,col+10) = [PrevDisp((t-1)*(nS+1)+11:t*(nS+1),j);PrevDisp((t-1)*(nS+1)+1:(t-1)*(nS+1)+9,j)]./UB(j);%Prev dispatch at t+9
%                 lockInputs((t-1)*nS+1:t*nS,col+11) = [PrevDisp((t-1)*(nS+1)+12:t*(nS+1),j);PrevDisp((t-1)*(nS+1)+1:(t-1)*(nS+1)+10,j)]./UB(j);%Prev dispatch at t+10
%                 lockInputs((t-1)*nS+1:t*nS,col+12) = [PrevDisp((t-1)*(nS+1)+13:t*(nS+1),j);PrevDisp((t-1)*(nS+1)+1:(t-1)*(nS+1)+11,j)]./UB(j);%Prev dispatch at t+11
%                 lockInputs((t-1)*nS+1:t*nS,col+13) = [PrevDisp((t-1)*(nS+1)+14:t*(nS+1),j);PrevDisp((t-1)*(nS+1)+1:(t-1)*(nS+1)+12,j)]./UB(j);%Prev dispatch at t+12
%                 lockInputs((t-1)*nS+1:t*nS,col+14) = [PrevDisp((t-1)*(nS+1)+15:t*(nS+1),j);PrevDisp((t-1)*(nS+1)+1:(t-1)*(nS+1)+13,j)]./UB(j);%Prev dispatch at t+13
%                 lockInputs((t-1)*nS+1:t*nS,col+15) = [PrevDisp((t-1)*(nS+1)+16:t*(nS+1),j);PrevDisp((t-1)*(nS+1)+1:(t-1)*(nS+1)+14,j)]./UB(j);%Prev dispatch at t+14
%                 lockInputs((t-1)*nS+1:t*nS,col+16) = [PrevDisp((t-1)*(nS+1)+17:t*(nS+1),j);PrevDisp((t-1)*(nS+1)+1:(t-1)*(nS+1)+15,j)]./UB(j);%Prev dispatch at t+15
%                 lockInputs((t-1)*nS+1:t*nS,col+17) = [PrevDisp((t-1)*(nS+1)+18:t*(nS+1),j);PrevDisp((t-1)*(nS+1)+1:(t-1)*(nS+1)+16,j)]./UB(j);%Prev dispatch at t+16
%                 lockInputs((t-1)*nS+1:t*nS,col+18) = [PrevDisp((t-1)*(nS+1)+19:t*(nS+1),j);PrevDisp((t-1)*(nS+1)+1:(t-1)*(nS+1)+17,j)]./UB(j);%Prev dispatch at t+17
%                 lockInputs((t-1)*nS+1:t*nS,col+19) = [PrevDisp((t-1)*(nS+1)+20:t*(nS+1),j);PrevDisp((t-1)*(nS+1)+1:(t-1)*(nS+1)+18,j)]./UB(j);%Prev dispatch at t+18
%                 lockInputs((t-1)*nS+1:t*nS,col+20) = [PrevDisp((t-1)*(nS+1)+21:t*(nS+1),j);PrevDisp((t-1)*(nS+1)+1:(t-1)*(nS+1)+19,j)]./UB(j);%Prev dispatch at t+19
%                 lockInputs((t-1)*nS+1:t*nS,col+21) = [PrevDisp((t-1)*(nS+1)+22:t*(nS+1),j);PrevDisp((t-1)*(nS+1)+1:(t-1)*(nS+1)+20,j)]./UB(j);%Prev dispatch at t+20
%                 lockInputs((t-1)*nS+1:t*nS,col+22) = [PrevDisp((t-1)*(nS+1)+23:t*(nS+1),j);PrevDisp((t-1)*(nS+1)+1:(t-1)*(nS+1)+21,j)]./UB(j);%Prev dispatch at t+21
%                 lockInputs((t-1)*nS+1:t*nS,col+23) = [PrevDisp((t-1)*(nS+1)+24:t*(nS+1),j);PrevDisp((t-1)*(nS+1)+1:(t-1)*(nS+1)+22,j)]./UB(j);%Prev dispatch at t+22
%                 lockInputs((t-1)*nS+1:t*nS,col+24) = [PrevDisp(t*(nS+1),j);PrevDisp((t-1)*(nS+1)+1:(t-1)*(nS+1)+23,j)]./UB(j);%Prev dispatch at t+23               
        end
    end
    col = col+3;
end

%prevent Nans
lockInputs(isnan(lockInputs)) = 0;

%convert from 0 to 1 scale to -1 to 1 scale
lockInputs = lockInputs*2-1;

if training %train it the first time through, use it later
    %if training separate into a training and testing set
    %use every ninth entry for testing so that you don't alias, but you
    %still have large enough training and testing sets
    
    %only include successful training sets 
%     lockInputs = lockInputs(1:8712*nS,:);
%     OnOffT = OnOffT(1:8712*nS,:);
    testi = [1:9:length(OnOffT(:,1))];
    %create vector of everything that isn't pulled out for testing
    traini = zeros(length(OnOffT(:,1))-length(testi),1);
    j =1;
    for i = 1:1:length(OnOffT(:,1))
        if rem(i,9) ~= 0
            traini(j) = i;
            j = j+1;
        end
    end
    lockInputs_train = lockInputs(traini,:);
    OnOff_train = OnOffT(traini,:);
    lockInputs_test = lockInputs(testi,:);
    OnOff_test = OnOffT(testi,:);
    
    %random order so that order doesn't effect output
    newOrder_train = randperm(length(traini));
    lockInputs_train = lockInputs_train(newOrder_train,:);
    OnOff_train = OnOff_train(newOrder_train,:);
    
    newOrder_test = randperm(length(testi));
    lockInputs_test = lockInputs_test(newOrder_test,:);
    OnOff_test = OnOff_test(newOrder_test,:);
    
    %create instance of class Neural_Network
    lockedNet = Neural_Network(length(lockInputs(1,:)),ngens,'classify',1);
    %train locked network
    disp('training ANN')
    tic
    [lockedNet, lockederror] = trainSingleLayer(lockedNet, OnOff_train,lockInputs_train);
    toc
    %test for accuracy
    Lockedtest = forward(lockedNet, lockInputs_test);
    Lockedtest(Lockedtest>=0.5) = true;
    Lockedtest(Lockedtest<0.5) = false;
    %accuracy is number correct of test set out of total
    acc = sum(Lockedtest==OnOff_test)./length(OnOff_test(:,1));
end
%use the network to find which generators to lock off or on
LockedNet = forward(lockedNet, lockInputs(1:nS,:))';
LockedNet(LockedNet>=0.5) = true;
LockedNet(LockedNet<0.5) = false;
%put the LockedNet, which is only for dispatchable components, back into
%the full Locked format
Locked = true(nS+1,nG);
Locked(2:end,gens) = LockedNet';

