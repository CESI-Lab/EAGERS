function  ThresholdSet(GenDisp,dt)
global Plant Threshold  DateSim 
 %% find threshold for any generators on/off at next step
 %% this is run in the Dipatchloop
IC = GenDisp(1,:);
EC = GenDisp(2,:);
Threshold =[];
QPall = Plant.Threshold.QP;
Organize = Plant.Threshold.Organize;

Outs = fieldnames(QPall);
[n,nG] = size(QPall.(Outs{1}).organize);
OnOff = true(1,nG);
for i = 1:1:nG
    OnOff(i) = Plant.Generator(i).Status;
end
n = n-1; %remove IC
dt = dt/n;
Time = linspace(1,n,n)*dt;
TurnOff = zeros(1,nG)+n;
TurnOn = zeros(1,nG)+n;
scaleCost = update_cost(Time/24+DateSim); %% All feedstock costs were assumed to be 1 when building matrices 
a = linspace(1,0,n+1)';
PredictDispatch = a*IC + (1-a)*EC;
% OriginalOptions = Plant.optimoptions;%in order for marginal cost to calculate correctly, you need to change the way time is done, but then change it back for everything else
% Plant.optimoptions.Horizon = Time(end);
% Plant.optimoptions.Resolution = dt;
% Plant.optimoptions.tspacing = 'constant';
marginCost = update_mc(Plant.Generator,PredictDispatch,scaleCost,Time);%the dispatch is whatever has been dispatched so far, except for the initial condition.
% Plant.optimoptions = OriginalOptions;
[QPall,Forecast] = updateMatrices(QPall,Organize,IC,Time,scaleCost,marginCost,EC);
%% Need to create range for forecast to be off by
OptDemands = Plant.optimoptions.Outputs;
for s = 1:1:length(OptDemands)
    r = n; % number of intervals in the uncertainty range to investigate
    a = linspace(.85,1.15,r);%the range of uncertainty --- scalers multiplied by the forecast
    Threshold.(OptDemands{s}).Range = Forecast.(OptDemands{s})'*a;
end
Threshold.Names = {};
for i = 1:1:nG
    Threshold.Names(i) = cellstr(Plant.Generator(i).Name);
end
include.E = {'CHP Generator', 'Electric Generator'};
include.C = {'Chiller'};
if ~Plant.optimoptions.sequential
    include.E = {'CHP Generator', 'Electric Generator', 'Chiller'};
end
Outs = fieldnames(QPall);
EnabledA = OnOff;
% EnabledA = IC>=LB;
EnabledB = EnabledA;
for s = 1:1:length(Outs)
    Threshold.(Outs{s}).Timestamp = DateSim+Time/24;
    Threshold.(Outs{s}).t = zeros(nG,1)+inf;
    Threshold.(Outs{s}).Upper = zeros(nG,1)+inf; % if there is a threshold that is unfeasible with current arrangement
    Threshold.(Outs{s}).Lower = zeros(nG,1)-inf; % if there is a threshold that is unfeasible with current arrangement

    gen = zeros(1,nG);
    dX_dt = zeros(1,nG);
    constCost =zeros(1,nG);
    QP = QPall.(Outs{s});
    Demand = Threshold.(Outs{s}).Range;
    if strcmp('E',Outs{s}) && ismember('H',OptDemands) %electricity & heat together
        DemandH = Threshold.H.Range;
    else DemandH = [];
    end
    for i = 1:1:nG
        if ismember(cellstr(Plant.Generator(i).Type),include.(Outs{s}))
            gen(i) = i;
            dX_dt(i) = Plant.Generator(i).OpMatB.Ramp.b(1); %ramp rate
            if isfield(Plant.Generator(i).OpMatB,'constCost')
                constCost(i) = Plant.Generator(i).OpMatB.constCost*dt;
            end
        end
    end
    gen = nonzeros(gen);
    dX_dt = dX_dt(gen);
    Threshold.(Outs{s}).Gen = gen;
    %% Find the next generators to turn on or off
    for j= 1:1:length(gen)
        i = gen(j);
        if OnOff(i) 
            r=[];
            if IC(i)>=LB(i) %only turn off if it has come all the way on
                r = find((GenDisp(2:end,i)<LB(i)),1,'first');%don't include the initial condition in the search
            else q = find((GenDisp(:,i)>=LB(i)),1,'first');
                if ~isempty(q)
                    r = q+find((GenDisp(q:end,i)<LB(i)),1,'first')-1;
                end
            end
            if ~isempty(r)
                if length(GenDisp(:,1))>=r+2
                    if GenDisp(r+2,i)<=LB(i) %it has to be off for more than one step to force it off
                        TurnOff(i) = r;
                    end
                end
            end
        else %generator is currently off
            r = find((GenDisp(2:end,i)>0),1,'first');%don't include the initial condition in the search
            if ~isempty(r)
                TurnOn(i) = r;
            end
        end
    end
    %% find lower threshold when things become infeasible and determine which generator to turn off
    % add up lower bounds of active generators (& subtract what can be dumped into energy storage)
    % if grid sellback then there is no lower threshold
    % pick generator with highest cost / next one to be shut down anyway
    
    %% find upper threshold when things become infeasible and determine which generator to turn on
    % add up upper bounds of active generators (& add what can be pulled from energy storage/grid)
    % pick generator with highest cost / next one to be shut down anyway
    
    %% find if any generators need to be turned on imediately due to ramping rates
    for j = 1:1:length(gen)
        RampDown = GenDisp(1,gen(j))-GenDisp(TurnOff(gen(j)),gen(j)); %amount it ramps down from now till shut off
        if RampDown/(dX_dt(j)*dt*TurnOff(gen(j)))> .95;
            Threshold.(Outs{s}).t(gen(j)) = 0; %turn off immediately
        else Threshold.(Outs{s}).t(gen(j)) = DateSim +(RampDown./dX_dt(j))/24; %time at which it needs to start shutting down to be off at next interval.
        end
        RampUp = GenDisp(TurnOn(gen(j)),gen(j))-GenDisp(1,gen(j));
        if RampUp/(dX_dt(j)*dt*TurnOn(gen(j)))> .95;
            Threshold.(Outs{s}).t(gen(j)) = 0; %turn on immediately
        else Threshold.(Outs{s}).t(gen(j)) = DateSim-(RampUp/dX_dt(j))/24;
        end
    end 
    %% find best time to change based on cost/value (create matrix of time and projected demand to interpolate later 
    tC1 = inf;
    tC2 = inf;
    ONtC = [];
    OFFtC = [];
    if  nnz(~isinf(TurnOn))>0
        [tC1,index] = sort(TurnOn(gen));
        ONtC = gen(index(tC1==tC1(1))); %list of generators turning on next
        tC1 = tC1+1;%add one to switch from number signifying next timestep, to index in GenDisp
    end
    if  nnz(~isinf(TurnOff))>0
        [tC2,index] = sort(TurnOff(gen));
        OFFtC = gen(index(tC2==tC2(1))); %list of generators turning off next
        tC2 = tC2+1;%add one to switch from number signifying next timestep, to index in GenDisp
    end
    Threshold.(Outs{s}).Off = [];
    Threshold.(Outs{s}).On = [];
    if min(tC1(1),tC2(1)) <= 3 %% generators are set to turn off/on in the next step
        if ~isempty(OFFtC) && tC2(1)<=3
            EnabledB(OFFtC) = 0;
            Threshold.(Outs{s}).Off = OFFtC;
        end
        if ~isempty(ONtC) && tC1(1)<=3
            EnabledB(ONtC) = 1;
            Threshold.(Outs{s}).On = ONtC;
        end
        [Threshold.(Outs{s}).Feasible,Threshold.(Outs{s}).Cost] = ThresholdFind(Outs{s},Demand,DemandH,QP,Organize,constCost,EnabledA,EnabledB);%figure out how to determine a threshold
    end
end
Threshold.Outputs = Outs;
Threshold.EnabledA = EnabledA;
Threshold.EnabledB = EnabledB;