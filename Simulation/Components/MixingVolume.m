function Out = MixingVolume(varargin) 
% a simple adiabatic mixing model mixing two or more inlet flows 
% Multiple (n+1) inlets: Flows (consisting of Temperature and flow rates of individual species) and P out 
% Three (3) outlets: Pressure, Flow , and Temperature 
% Many (n+2) states: Temperature, any species under consideration, and pressure
global Tags 
if length(varargin)==1 %first initialization
    block = varargin{1};
    block.Ru = 8.314472; % Universal gas constant in kJ/K*kmol
    block.Pdrop = 2; % presure drop across mixing volume (kPa)
    if ischar(block.SpeciesInit)
        block.SpeciesInit = ComponentProperty(block.SpeciesInit);
    end
    spec = fieldnames(block.SpeciesInit);
    spec = spec(~strcmp(spec,'T'));
    block.Scale = zeros(2+length(spec),1);
    block.IC = ones(length(block.Scale),1);
    block.Scale(1) = block.Tinit; 
    
    for i = 1:1:length(spec)
        if block.SpeciesInit.(spec{i}) ==0
            block.Scale(i+1,1) = NetFlow(block.SpeciesInit);
            block.IC(i+1,1) = 0;
        else
            block.Scale(i+1,1) = block.SpeciesInit.(spec{i});
        end
    end
    block.Scale(end) = 101 + block.Pdrop;
    block.UpperBound = inf*ones(length(block.Scale),1);
    block.LowerBound = zeros(length(block.Scale),1);
    block.spec = spec;
   
    block.InletPorts = {};
    for i = 1:1:block.inlets
        name = strcat('Inlet',num2str(i));
        block.(name) = [];
        block.InletPorts(end+1) = cellstr(name);
        block.(name).Saturation = [0,inf];
    end
    block.InletPorts(end+1) = {'Pout'};
    block.Pout.IC = 101;
    block.Pout.Saturation = [0,inf];
    block.Pout.Pstate = []; %identifies the state # of the pressure state if this block has one
    
    block.OutletPorts = {'Outlet';'Temperature';'Pin';};
    block.Outlet.IC = block.SpeciesInit;
    block.Outlet.IC.T = block.Tinit; 
    block.Temperature.IC = block.Outlet.IC.T;
    block.Pin.IC = block.Pout.IC+block.Pdrop;
    block.Pin.Pstate = length(block.Scale); %identifies the state # of the pressure state if this block has one
    
    block.P_Difference = {'Pin','Pout'};
    %no dMdP or mFlow (fixed pressure drop)
    Out = block;
elseif length(varargin)==2%% Have inlets connected, re-initialize
    block = varargin{1};
    Inlet =varargin{2};
    Inlet = checkSaturation(Inlet,block);
    n = block.inlets;     
    NetIn = {};
    for j = 1:1:length(block.spec)
        NetIn.(block.spec{j}) = 0;
    end
    inlets = fieldnames(Inlet);
    for i = 1:1:n
        spec2 = fieldnames(Inlet.(inlets{i}));
        spec2 = spec2(~strcmp('T',spec2));
        for j = 1:1:length(spec2)
            if ismember(spec2{j},block.spec)
                NetIn.(spec2{j}) = NetIn.(spec2{j}) + Inlet.(inlets{i}).(spec2{j});
            else
                block.spec(end+1) = spec2(j);
                NetIn.(spec2{j}) = Inlet.(inlets{i}).(spec2{j});
            end
        end
    end
    
    NetFlowIn = NetFlow(NetIn);
    block.Pfactor = NetFlowIn/block.Pdrop;
    
    H = zeros(n,1);
    Tin = zeros(n,1);
    for i = 1:1:n
        [~,H(i)] = enthalpy(Inlet.(inlets{i}));
        Tin(i) = Inlet.(inlets{i}).T;
    end
    
    Hin = sum(H);
    NetOut = NetIn;
    NetOut.T = mean(Tin);
    Terror = 1;
    while abs(Terror)>.01
        [~,Hout] = enthalpy(NetOut);
        Cp = SpecHeat(NetOut);
        Terror = (Hin-Hout)/(Cp*NetFlowIn);
        NetOut.T = NetOut.T+Terror;
    end
    block.Pout.IC = Inlet.Pout; 
    block.Pin.IC = Inlet.Pout+block.Pdrop;
    block.Outlet.IC = NetOut; 
    block.Temperature.IC = NetOut.T; 
    
    block.Scale = zeros(2+length(block.spec),1);
    block.IC = ones(length(block.Scale),1);
    block.Scale(1) = NetOut.T; 
    for i = 1:1:length(block.spec)
        if NetIn.(block.spec{i}) ==0
            block.Scale(i+1,1) = NetFlow(NetOut);
            block.IC(i+1,1) = 0;
        else
            block.Scale(i+1,1) = NetOut.(block.spec{i});
        end
    end
    block.Scale(end) = block.Pin.IC;
    block.UpperBound = inf*ones(length(block.Scale),1);
    block.LowerBound = zeros(length(block.Scale),1);
    Out = block;
else%running the model
    t = varargin{1};
    Y = varargin{2};
    Inlet = varargin{3};
    block = varargin{4};
    string1 = varargin{5};
    Inlet = checkSaturation(Inlet,block);
    Pin =Y(end);
    FlowOut = sum(max(0,Y(2:1+length(block.spec))));
    if FlowOut ==0
        scaleFlow = 1;
    else
        scaleFlow = (block.Pfactor*(Pin-Inlet.Pout))/FlowOut;%total flow out    
    end
    ActualOut.T = Y(1);
    for i = 1:1:length(block.spec)
        ActualOut.(block.spec{i}) = max(0,Y(i+1)*scaleFlow);
    end
    if strcmp(string1,'Outlet')
        Out.Pin = Pin;
        Out.Outlet = ActualOut;
        Out.Temperature = Out.Outlet.T;
        Tags.(block.name).MassFlow = MassFlow(ActualOut)/scaleFlow;
        Tags.(block.name).Temperature = Y(1);
    elseif strcmp(string1,'dY')
        dY = 0*Y;
        n = block.inlets;
        NetIn = {};
        for j = 1:1:length(block.spec)
            NetIn.(block.spec{j}) = 0;
        end
        inlets = fieldnames(Inlet);
        Hin = 0;
        for i = 1:1:n
            spec = fieldnames(Inlet.(inlets{i}));
            spec = spec(~strcmp('T',spec));
            for j = 1:1:length(spec)
                NetIn.(spec{j}) = NetIn.(spec{j}) + Inlet.(inlets{i}).(spec{j});
            end
            Hin = Hin + enthalpy(Inlet.(inlets{i}));
        end
        FlowIn = NetFlow(NetIn);
        FlowError = NetFlow(ActualOut) - FlowOut; 
        if FlowOut ==0
            Hout = 0;
            FlowOut = FlowIn;
        else
            scaleFlow2 = NetFlow(NetIn)/FlowOut ;
            Hout = enthalpy(ActualOut)*scaleFlow2/scaleFlow; %scale Hin and Hout to the same flow rate (no reactions so easy)
        end
        Cp = SpecHeat(ActualOut);
        dY(1) = (Hin-Hout)/(block.Vol*Cp*Pin).*Y(1)*block.Ru;
        %species
        for i = 1:1:length(block.spec)
%             dY(1+i) = (NetIn.(block.spec{i})-Y(i+1)).*Y(1)*block.Ru/(block.Vol*Pin);%change in stored mass of each species
            dY(1+i) = ((NetIn.(block.spec{i})/FlowIn-Y(i+1)/FlowOut)*FlowOut + FlowError*Y(i+1)/FlowOut)/(block.Vol*Pin).*Y(1)*block.Ru;%change in stored mass of each species
        end
        %pressure
        dY(end) = (FlowIn - NetFlow(ActualOut))*block.Ru*Y(1)/(block.Vol);
        Out = dY;
    end
end
end%Ends function MixingVolume