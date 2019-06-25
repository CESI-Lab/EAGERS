function [TurnOn, TurnOff] = ThresholdEnforce(Demand)
global OnOff GenAvailTime Threshold DateSim
TurnOff = [];
TurnOn = [];
S = Threshold.Outputs;
for s = 1:1:length(S)
    Thresh = Threshold.(S{s});
    Time = Thresh.Timestamp;
    %% make decision based on cost and forecasted demand
    if nnz(OnOff~=Threshold.EnabledB)>0 %still in A configuration, test if change is needed
        [n,r] = size(Thresh.Range); % # of timesteps (n) and number of demand steps (r)
        Time(end) = Time(end)+.0000001;%fix rounding error
        A = nnz(DateSim>Time)+1;
        if ~isinf(Thresh.Range(A,1))%make sure your range isn't infinite
            if Demand.(S{s})>Thresh.Range(A,end)
                b = r;
                disp('warning: Demand above +15% in thresholdEnforce')
            elseif Demand.(S{s})<Thresh.Range(A,1)
                b = 1;
                disp('warning: Demand below -15% in thresholdEnforce')
            else b = interp1(Thresh.Range(A,:),linspace(1,r,r),Demand.(S{s}));
            end
            if isnan(b)
                disp('thresholdEnforce')
            end
            CostFit = (ceil(b)-b)*Thresh.Cost(:,floor(b))+(b-floor(b))*Thresh.Cost(:,ceil(b));
            CostFit(isnan(CostFit)) = 0;
            [~,I] = min(CostFit);

            if I<=A && I<length(CostFit) %past optimal point to switch, and optimal is not at end
                for k = 1:1:length(Thresh.Off)
                    if OnOff(Thresh.Off(k))==1
                        TurnOff(end+1) = Thresh.Off(k);
                    end
                end
                for k = 1:1:length(Thresh.On)
                    if DateSim>=GenAvailTime(Thresh.On(k)) && OnOff(Thresh.On(k))==0
                        TurnOn(end+1) = Thresh.On(k);
                    end
                end
            end
        end
    end
    
    %% turn on or off generators by a certain deadline
    if DateSim>Thresh.t 
        
    end
    
    %% turn off generators to keep feasible
    if Demand.(S{s})<Thresh.Lower 
        %% should imediately trigger new dispatch
    end
    
    %% turn on generators to keep feasible
    if Demand.(S{s})>Thresh.Upper 
        %% should imediately trigger new dispatch
    end
    
end

TurnOn = unique(TurnOn);
TurnOff = unique(TurnOff);
D = datevec(DateSim);
if D(5)<10
    minute = strcat('0',num2str(D(5)));
else minute = num2str(D(5));
end
for k = 1:1:length(TurnOn)
    if OnOff(TurnOn(k))==1 %already on
        TurnOn = nonzeros(TurnOn(TurnOn~=TurnOn(k)));
    else
        name = Threshold.Names{TurnOn(k)};
        disp(['Turn On ',name,' at  ',num2str(D(4)),':',minute])
    end
end
for k = 1:1:length(TurnOff)
    if OnOff(TurnOff(k))==0 %already off
        TurnOff = nonzeros(TurnOff(TurnOff~=TurnOff(k)));
    else
        name = Threshold.Names{TurnOff(k)};
        disp(['Turn Off ',name,' at  ',num2str(D(4)),':',minute])
    end
end