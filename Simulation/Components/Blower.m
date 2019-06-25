function Out = Blower(varargin)
% A simple blower model that uses a sompressor map to relate PR, mass flow and efficiency
% Five (5) inlets: Inlet temperature, inlet species, outlet pressure, inlet pressure, power supply
% Two (2) outlets: Flow and speed (RPM)
% One (1) state: Speed
global Tags
if length(varargin)==1 % first initialization
    block = varargin{1};
    block.Ru = 8.314472; % Universal gas constant in kJ/K*kmol
    block.Mass = 1; % mass in kg
	block.InnerRadius = 0.05;
    block.OuterRadius = 0.15;
    block.Moment_Inertia = block.Mass/2*(block.InnerRadius^2 + block.InnerRadius^2);%Moment of Inertia for thick walled cylinder
    
    block.Scale = block.RPMdesign*2*pi/60; %speed in RPM converted to rad/s
    block.IC = ones(length(block.Scale),1);
    block.UpperBound = 1.5*block.Scale;
    block.LowerBound = .33*block.Scale;
    
    Dir=strrep(which('InitializeBlower.m'),fullfile('Components','Initialization','InitializeBlower.m'),'CompressorMaps');
    load(fullfile(Dir,block.Map));
    f = fieldnames(map);
    for i = 1:1:length(f)
        block.(f{i}) = map.(f{i});
    end
    Eff = map.Efficiency;
    Eff(Eff<0)=nan;
    block.minEfficiency = min(min(Eff));
    block.Efficiency(block.Efficiency<0) = block.minEfficiency;
    
    block.spec = fieldnames(block.InitialComposition);
    Mmass = MassFlow(block.InitialComposition);%kg/kmol
    Flow.T = block.Tdesign;
    for i = 1:1:length(block.spec)
        Flow.(block.spec{i}) = block.InitialComposition.(block.spec{i})*block.FlowDesign/Mmass;
    end
    
%     %ideal gas calculations
%     Cp = SpecHeat(Flow);
%     H1 = enthalpy(Flow);
%     T2s = Flow.T*(block.Pdesign)^((1.4 -1)/1.4);
%     Flow.T = T2s;
%     Cp = (SpecHeat(Flow)+Cp)/2;
%     Gamma = Cp/(Cp - block.Ru);
%     T2s = Flow.T*(block.Pdesign)^((Gamma -1)/Gamma);
%     Flow.T = T2s;
%     H2s = enthalpy(Flow);  
%     H2a = H1+(H2s-H1)/block.PeakEfficiency;
%     block.NominalPower = H2a - H1;
    
    %% engineering calculation (https://www.pdblowers.com/admin/uploads/rotary_positive_displacement_blowr_calculations_1.pdf)
    cfm = NetFlow(Flow)*1000*22.4*0.0353147*60; %convert kmol/s * 1000 mol/kmol * 22.4L/mol *0.0353147 ft^3/L * 60s/min
    psi = (block.Pdesign-1)*14.7;
    BHP = cfm*psi*.00436; 
    block.NominalPower = BHP/block.PeakEfficiency*0.7457;%converting horsepower to kW
    H1 = enthalpy(Flow);
    Cp = SpecHeat(Flow);
    Gamma = Cp/(Cp - block.Ru);
    Flow.T = Flow.T*(block.Pdesign)^((Gamma -1)/Gamma);
    error = 1;
    while abs(error)>1e-4
        H2 = enthalpy(Flow);
        error = (H1+block.NominalPower - H2)/(Cp*NetFlow(Flow));
        Flow.T = Flow.T + error;
    end

    %% set up dM/dP*Pout - C*Pin = mdot
    i = find(block.RPM>=1,1);
    j = ceil(length(block.NflowGMap(i,:))/2);
    P1 = block.PressRatio(i,j);
    P2 = block.PressRatio(i,j-1);
    M1 = block.NflowGMap(i,j);
    M2 = block.NflowGMap(i,j-1);
    block.mFlow = block.FlowDesign;
    dMdP = block.FlowDesign*(M2 - M1)/(101*(block.Pdesign-1)*(P2 - P1));
    C = (dMdP*101*block.Pdesign-block.mFlow)/101;
    block.dMdP = [dMdP, C]; 

    %%
    block.InletPorts = {'Temperature','Species','Pin','Pout','Power'};
    block.Temperature.IC = block.Tdesign; 
    block.Temperature.Saturation = [0,inf];
    block.Species.IC = block.InitialComposition; 
    block.Species.Saturation = [0,1];
    block.Pin.IC = 101;%in kPa
    block.Pin.Saturation = [0,inf];
    block.Pin.Pstate = []; %identifies the state # of the pressure state if this block has one
    block.Pout.IC = block.Pin.IC*block.Pdesign;%in kPa
    block.Pout.Saturation = [0,inf];
    block.Pout.Pstate = []; %identifies the state # of the pressure state if this block has one
    block.Power.IC = block.NominalPower;
    block.Power.Saturation = [0,inf];
    
    block.OutletPorts = {'Outlet','Speed'};
    block.Outlet.IC = Flow;
    block.Speed.IC = block.RPMdesign;
    
    block.P_Difference = {'Pout','Pin'};
    
    Out = block;
elseif length(varargin)==2 %% Have inlets connected, re-initialize
    block = varargin{1};
    Inlet = varargin{2};
    block.spec = fieldnames(Inlet.Species);
    Inlet = checkSaturation(Inlet,block);
    %problem is to find speed and flow that matches the PR and power specified by input
    RPM = block.Scale(1)*60/(2*pi); %convert rad/s back to RPM
    [Power,Flow] = OpPoint(RPM,Inlet,block);
    block.NominalPower = Power;

    FlowIn = Flow;
    FlowIn.T = Inlet.Temperature;
    H2a = enthalpy(FlowIn)+block.NominalPower;
    Cp = SpecHeat(Flow);
    tC = (Cp*NetFlow(Flow));
    errorT =1;
    while abs(errorT)>.01 %solve for the correct outlet temperature
        Hout = enthalpy(Flow);
        errorT = (H2a - Hout)/tC;
        Flow.T = Flow.T + errorT;
    end
    
    block.IC = (RPM*2*pi/60)/block.Scale; %speed in RPM converted to rad/s
    block.Outlet.IC = Flow;
    block.Speed.IC = RPM;
    block.mFlow = MassFlow(Flow);
    
    %% set up dM/dP*Pout - C*Pin = mdot
    Mscale = block.FlowDesign*(Inlet.Pin/block.P0map)/(MassFlow(Inlet.Species)*(Inlet.Temperature/block.Tdesign)^.5);%scalar maping for mass flow
    CompPR = (Inlet.Pout/Inlet.Pin -1)/(block.Pdesign - 1) + 1;%normalized compressor pressure ratio
    nRPM = RPM/(block.RPMdesign*(Inlet.Temperature/block.Tdesign)^.5);%normalized RPM
    i2 = find(block.RPM >=nRPM,1,'first');
    if isempty(i2)
        i2 = length(block.RPM);
    elseif i2==1
        i2 =2;
    end
    i1 = i2-1;
    r = (nRPM-block.RPM(i1))/(block.RPM(i2)-block.RPM(i1));
    j12 = find(block.PressRatio(i1,:)>=CompPR,1,'first');
    if isempty(j12)
        j12 = length(block.PressRatio(i1,:));
    elseif j12==1
        j12 =2;
    end
    j11 = j12-1;
    j22 = find(block.PressRatio(i2,:)>=CompPR,1,'first');
    if isempty(j22)
        j22 = length(block.PressRatio(i2,:));
    elseif j22==1
        j22 = 2;
    end
    j21 = j22-1;
    P1 = block.PressRatio(i1,j11)*r+ block.PressRatio(i2,j21)*(1-r);
    P2 = block.PressRatio(i1,j12)*r+ block.PressRatio(i2,j22)*(1-r);
    P1 = Inlet.Pin*((P1 -1)*(block.Pdesign -1) + 1);
    P2 = Inlet.Pin*((P2 -1)*(block.Pdesign -1) + 1);
    M1 = block.NflowGMap(i1,j11)*r+ block.NflowGMap(i2,j21)*(1-r);
    M2 = block.NflowGMap(i1,j12)*r+ block.NflowGMap(i2,j22)*(1-r);
    if M1 == M2
        if j12<length(block.NflowGMap(i1,:))
            M2 = block.NflowGMap(i1,j12+1)*r+ block.NflowGMap(i2,j22+1)*(1-r);
        else
            M1 = block.NflowGMap(i1,j11-1)*r+ block.NflowGMap(i2,j21-1)*(1-r);
        end
    end
    dMdP = Mscale*(M2 - M1)/(P2 - P1);
    C =  (block.dMdP(1,1)*Inlet.Pout-block.mFlow)/Inlet.Pin;
    block.dMdP = [dMdP, C]; 
    
    Tags.(block.name).RPM = RPM;
    Tags.(block.name).Flow = Flow;
    Tags.(block.name).Power = Inlet.Power;
    Tags.(block.name).MassFlow = MassFlow(Flow);
    Tags.(block.name).Temperature = Flow.T;
    Out = block;
else%running the model
    t = varargin{1};
    Y = varargin{2};
    Inlet = varargin{3};
    block = varargin{4};
    string1 = varargin{5};
    
    Inlet = checkSaturation(Inlet,block);
    RPM = Y(1)*60/(2*pi); %convert rad/s back to RPM
    [Power,Flow] = OpPoint(RPM,Inlet,block);
    NetPower = Inlet.Power - Power;
    if strcmp(string1,'Outlet')
        Out.Outlet = Flow;
        Out.Speed = RPM;
        Tags.(block.name).RPM = RPM;
        Tags.(block.name).Flow = Flow;
        Tags.(block.name).Power = Inlet.Power;
        Tags.(block.name).PR = Inlet.Pout/Inlet.Pin;
        Tags.(block.name).MassFlow = MassFlow(Flow);
        Tags.(block.name).Temperature = Flow.T;
        Tags.(block.name).nMflow = MassFlow(Flow)/block.FlowDesign;
        Tags.(block.name).Pressure = Inlet.Pout;
    elseif strcmp(string1,'dY')
        dY = 0*Y;
        dY(1) = NetPower/(block.Moment_Inertia*Y(1)); %dw/dt = P/(w*l) units of rad/s
        Out = dY;
    end
end
end %Ends function Blower


function [Power,Flow] = OpPoint(RPM,Inlet,block)
global Tags
nT = Inlet.Temperature/block.Tdesign;%square root of the normalized T
nRPM = min(block.RPM(end),max(block.RPM(1),RPM/(block.RPMdesign*nT^.5)));%normalized RPM
CompPR = (Inlet.Pout/Inlet.Pin -1)/(block.Pdesign - 1) + 1;%normalized compressor pressure ratio
Mmass = MassFlow(Inlet.Species);%kg/kmol
Mscale = block.FlowDesign*(Inlet.Pin/block.P0map)/(Mmass*(Inlet.Temperature/block.Tdesign)^.5);%scalar maping for mass flow
%% Find table indecies (i,j) surrounding actual value
% s and r correspond to the fractional coordinate position 
i2 = find(block.RPM >=nRPM,1,'first');
if isempty(i2)
    i2 = length(block.RPM);
elseif i2==1
    i2 =2;
end
i1 = i2-1;
s =(nRPM - block.RPM(i1))/(block.RPM(i2) - block.RPM(i1));
PRVec = block.PressRatio(i1,:)*(1-s) + block.PressRatio(i2,:)*s;
if CompPR>max(PRVec) && nRPM<1
    CompPR = (Inlet.Pout/Inlet.Pin -1)*nRPM/(block.Pdesign - 1) + 1;%renormalized compressor pressure ratio
end
Beta = interp1(PRVec,block.Beta,CompPR,'spline');
Nflow = interp2(block.Beta,block.RPM,block.NflowGMap,Beta,nRPM,'spline');
Eff = interp2(block.Beta,block.RPM,block.Efficiency,Beta,nRPM,'spline')*block.PeakEfficiency;
NetFlowOut = Nflow*Mscale;

Flow.T = Inlet.Temperature;
for i = 1:1:length(block.spec)
    Flow.(block.spec{i}) = Inlet.Species.(block.spec{i})*NetFlowOut;
end

% %Ideal gas calculation
% Cp = SpecHeat(Flow);
% H1 = enthalpy(Flow);
% T2s = Flow.T*(Inlet.Pout/Inlet.Pin)^((1.4 -1)/1.4);
% Flow.T = T2s;
% Cp = (SpecHeat(Flow)+Cp)/2;
% Gamma = Cp/(Cp - block.Ru);
% T2s = Flow.T*(Inlet.Pout/Inlet.Pin)^((Gamma -1)/Gamma);
% Flow.T = T2s;
% H2s = enthalpy(Flow);  
% H2a = H1+(H2s-H1)/Eff;
% Power = (H2a - H1);

%% engineering calculation (https://www.pdblowers.com/admin/uploads/rotary_positive_displacement_blowr_calculations_1.pdf)
cfm = NetFlow(Flow)*1000*22.4*0.0353147*60; %convert kmol/s * 1000 mol/kmol * 22.4L/mol *0.0353147 ft^3/L * 60s/min
psi = (Inlet.Pout - Inlet.Pin)*14.7/101.325;
BHP = cfm*psi*.00436; 
Power = BHP/Eff*0.7457;%converting horsepower to kW
H1 = enthalpy(Flow);
Cp = SpecHeat(Flow);
Gamma = Cp/(Cp - block.Ru);
Flow.T = Flow.T*(Inlet.Pout/Inlet.Pin)^((Gamma -1)/Gamma);
error = 1;
while abs(error)>1e-4
    H2 = enthalpy(Flow);
    error = (H1+Power - H2)/(Cp*NetFlow(Flow));
    Flow.T = Flow.T + error;
end

Tags.(block.name).NRPM = nRPM;
Tags.(block.name).Beta = Beta;
Tags.(block.name).Eff = Eff;
Tags.(block.name).Nflow = Nflow;
end %Ends function OpPoint