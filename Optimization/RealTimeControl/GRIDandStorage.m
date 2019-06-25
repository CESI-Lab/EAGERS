function Out = GRIDandStorage(Demand,setGen,stor,utility)
global targetCharge targetDischarge LBmpc UBmpc
setStor =setGen(stor);
netStor = sum(setStor);
nStor = length(stor);
setGrid =setGen(utility);
netGrid  = sum(setGrid);
nGrid = length(setGrid);

if nGrid==0 && nStor>0 %no grid, all error is absorbed by storage
    i=1;
    newDischarge = max(Demand-netStor,0);
    newCharge = max(netStor-Demand,0);
    while newDischarge>1e-8 && i<=nStor%distrubute discharge amongst energy storage
        a = min(newDischarge,UBmpc(stor(i))-setStor(i));
        setStor(i) = setStor(i)+a;
        newDischarge = newDischarge-a;
        i = i+1;
    end
    while newCharge>1e-8 && i<=nStor %distrubute charging amongst energy storage
        a = min(newCharge,-LBmpc(stor(i))+setStor(i));
        setStor(i) = setStor(i)-a;
        newCharge = newCharge-a;
        i = i+1;
    end
elseif isempty(nStor) && ~isempty(nGrid) %no storage, grid handles error
    netGrid = Demand;% split grid to buying and selling if necessary
else
    %split the demand which must be met instantaneously by grid or energy
    %storage between grid and energy storage

    %increase discharge towards target if this reduces grid purchases
    if netGrid>1e-4 && netStor>0 && netStor<sum(targetDischarge(1:nStor))% if buying power, discharging, and using less than the average expected discharge
        newDischarge = min(netGrid,sum(targetDischarge(1:nStor))-netStor); %increase discharge until a) grid purchase =0 or b) discharge = target
        i = 1;
        while newDischarge>0 && i<=nStor%distrubute discharge amongst energy storage
            a = min(newDischarge,targetDischarge(i));
            setStor(i) = setStor(i)+a;
            newDischarge = newDischarge-a;
            i = i+1;
        end
        netStor = sum(setStor);
        netGrid= Demand-netStor; %reduce purchased power
    end

    %decrease charge towards target & beyond if this reduces grid purchases
    if netGrid>1e-4 && netStor<0 && -netStor>sum(targetCharge(1:nStor))% if buying power, charging, and charging more than expected
        lessCharge = max(netGrid,-netStor-sum(targetCharge(1:nStor))); %decrease charging until a) grid purchase =0 or b) charge = target
        i = 1;
        while lessCharge>0 && i<=nStor%distrubute discharge amongst energy storage
            a = min(lessCharge,-setStor(i)-targetCharge(i));
            setStor(i)= setStor(i)+a;
            lessCharge = lessCharge-a;
            i = i+1;
        end
        netStor = sum(setStor);
        netGrid= Demand-netStor; %reduce purchased power
        i=nStor;
        while netGrid>sum(setGrid) && netStor<0 && i>0%buying more than planned, reduce charging to zero
            a = min(netGrid-sum(setGrid),-setStor(i));
            setStor(i) = setStor(i)+a;
            netGrid = netGrid-a;
            i=i-1;
        end
        netStor = sum(setStor);
    end
    error = Demand-netGrid-netStor;
    if error>1e-4  && netGrid<0 %less generation/more demand than expected & increasing discharge towards target was not enough
        if -netGrid>error%reduce amount of power sold
            netGrid = netGrid+error;
        else %set power sold to zero
            netGrid =0;
        end
        error = Demand-netGrid-netStor;
    end

    if error>1e-4 && netGrid>=0 %reducing sold power was not enough
        if sum(UBmpc(stor))<(Demand-netGrid) %maximum discharge
            netGrid = Demand-sum(UBmpc(stor));
            setStor = UBmpc(stor);
        else %buy/sell the same & discharge more
            newDischarge = Demand-netGrid-netStor;
            i = 1;
            while newDischarge>0 && i<=nStor%distrubute additional discharge amongst energy storage
                a = min(newDischarge,UBmpc(stor(i))-setStor(i));
                setStor(i) = setStor(i)+a;
                newDischarge = newDischarge-a;
                i = i+1;
            end
        end
        netStor = sum(setStor);
        error = Demand-netGrid-netStor;
    end
    if error<-1e-4 %more generation/less demand than expected
        if netStor>sum(targetCharge(1:nStor))
            sFrac = setStor/netStor;
            setStor = setStor + max(error*sFrac, targetCharge(1:nStor)-setStor);
            netStor = sum(setStor);
            error = Demand-netGrid-netStor;
        end
    end
    if error<-1e-4 %more generation/less demand than expected
        if netGrid>0 && (Demand-netStor)>=0%currently buying power & buying less balances equation
            %buy less power (charge/discharge the same)
        elseif netGrid>0 %buying power, but prefer to buy less & charge less
            newCharge = min(-(Demand-netStor),(-sum(LBmpc(stor))+netStor));
            i = 1;
            while newCharge>0 && i<=nStor%distrubute discharge amongst energy storage
                a = min(newCharge,(-LBmpc(stor(i))+setStor(i)));
                setStor(i) = setStor(i)-a;
                newCharge = newCharge-a;
                i = i+1;
            end
        elseif sum(LBmpc(stor))>Demand-netGrid %selling power, a) increase charging up to limit, b) increase charging up to limit && increase selling but exceed maximum Charge rate
            setStor = LBmpc(stor);
        end
        netStor = sum(setStor);
        netGrid = Demand-netStor;
    end
%         if netGrid<0
%             i=1;
%             while netGrid<0 && i<=nStor
%                 a=min(-netGrid,(-LBmpc(stor(i))+setStor(i)));
%                 setStor(i) = setStor(i)-a;%charge storage instead of exporting power
%                 netGrid = netGrid+a;
%                 i=i+1;
%             end
%         end  
end
if abs(sum(setGrid))>0
    setGrid = netGrid*setGrid/abs(sum(setGrid));%this makes sure that the change is made to the grid component that is being used.
elseif ~isempty(utility)
    setGrid(1) = netGrid; %may run into issues with multiple grids
end
Out = [setStor, setGrid];