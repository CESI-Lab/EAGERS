function Out = Oxidizer(varargin)
% an oxidizer with n inlet flows 
% Complete combustion is assumed for all fuels (if enough O2 present)
% two or more (2+) inlets: Flows and outletl pressure
% Two (2) outlets: Flow out and pressure at the inlet
% Many (2+n) states: Outlet Temperature, Species, Inlet Pressure
global Tags
if length(varargin)==1 % first initialization
    block = varargin{1};
    block.Ru = 8.314472; % Universal gas constant in kJ/K*kmol
    % CH4+1.5O2 --> CO+ 2H2O
    % CO + .5O2 --> CO2
	% H2 + .5O2 --> H2O
    
    block.spec = {'CH4';'CO';'CO2';'H2';'H2O';'O2'};
    block.Scale = block.InitialFlowOut.T;
    SpecNew = fieldnames(block.InitialFlowOut);
    SpecAll = unique([block.spec;SpecNew]);
    block.spec = SpecAll(~strcmp('T',SpecAll));
    for i = 1:1:length(block.spec)
        if ~ismember(block.spec{i},SpecNew)
            block.InitialFlowOut.(block.spec{i}) = 0;
        end
        block.Scale(end+1,1) = block.InitialFlowOut.(block.spec{i});
    end
    
    block.Pdrop = 2; % presure drop across mixing volume (kPa)
    block.Scale(end+1,1) = 101+block.Pdrop;
    block.IC = ones(length(block.Scale),1);
    block.UpperBound = inf*ones(length(block.Scale),1);
    block.LowerBound = zeros(length(block.Scale),1);
    
    block.InletPorts = {};
    for j = 1:1:block.inlets
        name = strcat('Inlet',num2str(j));
        block.InletPorts(end+1) = cellstr(name);
        block.(name).IC.T = 800;
        block.(name).Saturation = [0,inf];
    end
    block.InletPorts(end+1) = {'Pout'};
    
    block.Pout.type = 'in';
    block.Pout.IC = 101;
    block.Pout.Saturation = [0,inf];
    block.Pout.Pstate = []; %identifies the state # of the pressure state if this block has one
    
    block.OutletPorts = {'Flow';'Pin';'MeasureT'};
    block.Flow.IC = block.InitialFlowOut;
    block.Pin.IC = block.Pout.IC+block.Pdrop;
    block.Pin.Pstate = length(block.Scale); %identifies the state # of the pressure state if this block has one
    block.MeasureT.IC = block.InitialFlowOut.T;
    
    block.P_Difference = {'Pin','Pout'};
    %no dMdP or mFlow (fixed pressure drop)
    Out = block;
elseif length(varargin)==2 %% Have inlets connected, re-initialize
    block = varargin{1};
    Inlet = varargin{2};
    Inlet = checkSaturation(Inlet,block);
    n = block.inlets;
    inlets = fieldnames(Inlet);
    Pin = Inlet.Pout+block.Pdrop;
    %mix flows
    NetIn =[];
    H_in = 0;
    for j = 1:1:n
        H_in = H_in + enthalpy(Inlet.(inlets{j}));
        SpecNew = fieldnames(Inlet.(inlets{j}));
        SpecAll = unique([block.spec;SpecNew]);
        block.spec = SpecAll(~strcmp('T',SpecAll));
        for i = 1:1:length(block.spec)
            if ~isfield(NetIn,block.spec{i})
                NetIn.(block.spec{i}) = 0;
            end
            if ismember(block.spec{i},SpecNew)
                NetIn.(block.spec{i}) = NetIn.(block.spec{i}) + Inlet.(inlets{j}).(block.spec{i});
            end
        end
    end
    
    ReactMix.T = block.Scale(1);
    R.CH4 = NetIn.CH4;
    R.H2 = NetIn.H2;
    if ReactMix.T>1700
        Kp = 0.4589*(ReactMix.T/3000)^14.772;
    else
        Kp = max(0,(ReactMix.T-300)*(6.78e-5/1400));
    end
    X_CO_X_CO2 = Kp/((Pin/101)^.5*(NetIn.O2/NetFlow(NetIn))); %estimate ratio of CO to CO2
    R.CO = ((NetIn.CO+R.CH4)-X_CO_X_CO2*NetIn.CO2)/(1+X_CO_X_CO2);

    sumR = 1.5*R.CH4 + 0.5*R.CO + 0.5*R.H2;
    phi = sumR/NetIn.O2;% O2 necessary/O2 supplied
    r = fieldnames(R);
    lean = min(1,NetIn.O2/sumR);
    for j = 1:1:length(r)
        R.(r{j}) = R.(r{j})*lean; %rich combustion limited by O2
    end
    
    for i = 1:1:length(block.spec)
        if strcmp(block.spec{i},'CH4')
            ReactMix.CH4 = NetIn.CH4 - R.CH4;
        elseif strcmp(block.spec{i},'CO')
            ReactMix.CO = NetIn.CO + R.CH4- R.CO;
        elseif strcmp(block.spec{i},'CO2')
            ReactMix.CO2 = NetIn.CO2 + R.CO;
        elseif strcmp(block.spec{i},'H2')
            ReactMix.H2 = NetIn.H2  - R.H2;
        elseif strcmp(block.spec{i},'H2O')
            ReactMix.H2O = NetIn.H2O + 2*R.CH4 + R.H2;
        elseif strcmp(block.spec{i},'O2')
            ReactMix.O2 = NetIn.O2 - 1.5*R.CH4 - .5*R.CO - .5*R.H2;
        else
            ReactMix.(block.spec{i}) = NetIn.(block.spec{i});
        end
    end
    Flow = NetFlow(ReactMix);
    Terror = 100;
    while abs(Terror)>1e-2
        Cp = SpecHeat(ReactMix);
        Terror = (H_in-enthalpy(ReactMix))/(Cp*Flow);
        ReactMix.T = ReactMix.T+Terror;
    end
    block.Flow.IC = ReactMix;
    block.Pfactor = NetFlow(ReactMix)/block.Pdrop;
    %no dMdP or mFlow to update (fixed pressure drop)
    block.Pout.IC = Inlet.Pout; 
    block.Pin.IC = Inlet.Pout+block.Pdrop;
    block.MeasureT.IC = ReactMix.T;
    block.Scale = ReactMix.T;
    block.IC = 1;
    for i = 1:1:length(block.spec)
        if ReactMix.(block.spec{i})==0
            block.IC(end+1,1) = ReactMix.(block.spec{i})/Flow;
            block.Scale(end+1,1) = Flow;
        else
            block.IC(end+1,1) = 1;
            block.Scale(end+1,1) = ReactMix.(block.spec{i});
        end
    end
    block.Scale(end+1,1) = Inlet.Pout+block.Pdrop;
    block.IC(end+1,1) = 1;
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
    Pin = Y(end);
    %merge flows
    Hin = 0;
    NetIn =[];
    n = block.inlets;
    inlets = fieldnames(Inlet);
    for j = 1:1:n
        Hin = Hin + enthalpy(Inlet.(inlets{j}));
        Spec = fieldnames(Inlet.(inlets{j}));
        for i = 1:1:length(block.spec)
            if ~isfield(NetIn,block.spec{i})
                NetIn.(block.spec{i}) = 0;
            end
            if ismember(block.spec{i},Spec)
                NetIn.(block.spec{i}) = NetIn.(block.spec{i}) + Inlet.(inlets{j}).(block.spec{i});
            end
        end
    end
    ReactMix.T = Y(1);
    %% 3 reaction:
    % CH4 + 1.5O2 --> CO + 2H2O
    % CO + .5O2 --> CO2
    % H2 + .5O2 --> H2O
    R.CH4 = NetIn.CH4;
    R.H2 = NetIn.H2;
    if ReactMix.T>1700
        Kp = 0.4589*(ReactMix.T/3000)^14.772;
    else
        Kp = max(0,(ReactMix.T-300)*(6.78e-5/1400));
    end
    X_CO_X_CO2 = Kp/((Pin/101)^.5*(NetIn.O2/NetFlow(NetIn))); %estimate ratio of CO to CO2
    R.CO = ((NetIn.CO+R.CH4)-X_CO_X_CO2*NetIn.CO2)/(1+X_CO_X_CO2);

    sumR = 1.5*R.CH4 + 0.5*R.CO + 0.5*R.H2;
    phi = sumR/NetIn.O2;

    scaleFlow = (block.Pfactor*(Pin-Inlet.Pout))/sum(max(0,Y(2:1+length(block.spec))));%total flow out
    if scaleFlow<=0
        scaleFlow=1;
    end
    ActualOut.T = Y(1);
    for i = 1:length(block.spec)
        ActualOut.(block.spec{i}) = max(0,Y(1+i))*scaleFlow;
    end

    if strcmp(string1,'Outlet')
        Out.Flow = ActualOut;
        Out.Pin = Pin;
        Out.MeasureT = ActualOut.T;
        Tags.(block.name).EquivelanceRatio = phi;
        Tags.(block.name).Temperature = Y(1);
        Tags.(block.name).MassFlow = MassFlow(Out.Flow);
    elseif strcmp(string1,'dY')
        dY = 0*Y;

        r = fieldnames(R);
        if phi>1 %rich combustion
            for i = 1:1:length(r)
                R.(r{i}) = R.(r{i})/phi; %rich combustion limited by O2
            end
        end
        for i = 1:1:length(block.spec)
            if strcmp(block.spec{i},'CH4')
                ReactMix.CH4 = NetIn.CH4 - R.CH4;
            elseif strcmp(block.spec{i},'CO')
                ReactMix.CO = NetIn.CO + R.CH4- R.CO;
            elseif strcmp(block.spec{i},'CO2')
                ReactMix.CO2 = NetIn.CO2 + R.CO;
            elseif strcmp(block.spec{i},'H2')
                ReactMix.H2 = NetIn.H2  - R.H2;
            elseif strcmp(block.spec{i},'H2O')
                ReactMix.H2O = NetIn.H2O + 2*R.CH4 + R.H2;
            elseif strcmp(block.spec{i},'O2')
                ReactMix.O2 = NetIn.O2 - 1.5*R.CH4 - .5*R.CO - .5*R.H2;
            else
                ReactMix.(block.spec{i}) = NetIn.(block.spec{i});
            end
        end
        ReactFlow = NetFlow(ReactMix);
        FlowOut = sum(max(0,Y(2:1+length(block.spec))));
        FlowError = NetFlow(ActualOut) - FlowOut;
        %temperature
        if FlowOut<=0
            scaleFlow2 = 1;
        else
            scaleFlow2 = NetFlow(ReactMix)/FlowOut;
        end
        Hout = enthalpy(ActualOut)*scaleFlow2/scaleFlow; %scale Hin and Hout to the same flow rate (no reactions so easy)
        Cp = SpecHeat(ActualOut);
        dY(1) = (Hin-Hout)/(block.Vol*Cp*Pin).*Y(1)*block.Ru;
        for i = 1:1:length(block.spec)
    %         dY(1+i) = (ReactMix.(block.spec{i})-Y(i+1))*block.Ru*Y(1)/(block.Vol*Pin);%change in stored mass of each species
            dY(1+i) = ((ReactMix.(block.spec{i})/ReactFlow-Y(i+1)/FlowOut)*FlowOut + FlowError*max(0,Y(i+1))/FlowOut)*block.Ru*Y(1)/(block.Vol*Pin);%change in stored mass of each species
        end
        dY(2+i) = (ReactFlow - NetFlow(ActualOut))*block.Ru*Y(1)/(block.Vol); %change in Pressure
        Out = dY;
    end
end
end%Ends function Oxidizer