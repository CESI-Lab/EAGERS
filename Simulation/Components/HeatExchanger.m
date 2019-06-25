function Out = HeatExchanger(varargin)
%Nodal Heat Exchanger model with 3 states per node: Temperature hot, Plate temperature, temperature cold
% Four (4) inlets: {'Cold Flow','Hot Flow','Cold Pout','Hot Pout'}
% Four (4) outlets: {'Cold Flow','Hot Flow','Cold Pin','Hot Pin'}
% Three* (3*n) states: Cold, solid, and hot flow temperature states for each node
global Tags
if length(varargin)==1 % first initialization
    block = varargin{1};
    block.Ru = 8.314472; % Universal gas constant in kJ/K*kmol
    %% Load mask parameters (flow direction)
    block.nodes = block.rows*block.columns;
    block = FlowDir(block,2);
    
    if ischar(block.ColdSpecIn)
        block.ColdSpecIn = ComponentProperty(block.ColdSpecIn);
    end
    block.spec1 = fieldnames(block.ColdSpecIn);
    block.spec1 = block.spec1(~strcmp(block.spec1,'T'));
    for i = 1:1:length(block.spec1)
        Inlet.Flow1.(block.spec1{i}) = block.ColdSpecIn.(block.spec1{i});
    end
    Inlet.Flow1.T = block.Cold_T_init;
    
    if ischar(block.HotSpecIn)
        block.HotSpecIn = ComponentProperty(block.HotSpecIn);
    end
    block.spec2 = fieldnames(block.HotSpecIn);
    block.spec2 = block.spec2(~strcmp(block.spec2,'T'));
    for i = 1:1:length(block.spec2)
        Inlet.Flow2.(block.spec2{i}) = block.HotSpecIn.(block.spec2{i});
    end
    Inlet.Flow2.T = block.Hot_T_init;
    
    %% Heat Exchanger
    block.PdropCold = 1;                                    % (kPa) pressure drop
    block.PdropHot = 1;                                    % (kPa) pressure drop
    block.h_conv= 50;                                       %Convective heat coefficient (W/m^2*K)
    block.t_Plate = .001;                                   %Thickness of Heat Exchanger plate (m)
    block.Length = block.Vol^(1/3);                         %Length of heat Exchanger (m)
    block.Width = block.Vol^(1/3);                                   %Length of heat Exchanger (m)
    block.Solid_SpecHeat=.48;                               %Specific Heat of solid material (kJ/kg*K)
    block.Solid_CondCoef = 0;%15;                              %Conduction coefficient (Steel) (W/m*K)
    block.Solid_Density=8055;                               %Density of solid material (kg/m^3)
    block.Vol_Solid = min(1/3*block.Vol, block.Mass/block.Solid_Density);
    block.Vol_Cold = (block.Vol-block.Vol_Solid)/2;
    block.Vol_Hot = (block.Vol-block.Vol_Solid)/2;
    
    Target = ComponentProperty(block.Target);
    [block.Area,~,~] = findArea(Inlet.Flow1,Inlet.Flow2,block.h_conv,Target,block.sizemethod);
    block.Convection = block.h_conv*block.Area/block.nodes/1000; %all heat transfer coefficients converted to kW/K: thus Q = C*(T1-T2) is in kW
    %horizontal
    block.L_node = block.Length/block.columns;
    block.W_node =block.Width/block.rows;
    block.AcondLR = block.Vol_Solid/block.Width/block.columns;
    block.AcondPN = block.Vol_Solid/block.Length/block.rows;
    block.ConductionPN = block.Solid_CondCoef*block.AcondPN/(block.L_node/2)/1000; %heat transfer coefficient between previous and next node of plate
    block.ConductionLR  = block.Solid_CondCoef*block.AcondLR/(block.W_node/2)/1000; %heat transfer coefficient between left and right adjacent nodes of oxidant plate
    
    mdot_Cp(1) = SpecHeat(Inlet.Flow1).*NetFlow(Inlet.Flow1)/block.rows;
    if strcmp(block.direction ,'crossflow')
        mdot_Cp(2) = SpecHeat(Inlet.Flow2).*NetFlow(Inlet.Flow2)/block.columns;
    else
        mdot_Cp(2) = SpecHeat(Inlet.Flow2).*NetFlow(Inlet.Flow2)/block.rows;
    end
    [block.Scale , block.HTcond, block.HTconv] = SteadyTemps(block,mdot_Cp,[Inlet.Flow1.T,Inlet.Flow2.T]);
    block.Scale(end+1:end+2,1) = [101+block.PdropCold;101+block.PdropHot;];
    block.IC = ones(3*block.nodes+2,1);
    block.UpperBound = inf*ones(3*block.nodes+2,1);
    block.LowerBound = zeros(3*block.nodes+2,1);
    
    Flow1 = Inlet.Flow1;
    Flow1.T =  mean(block.Scale(block.Flow1Dir(:,end),1));
    Flow2 = Inlet.Flow2;
    Flow2.T = mean(block.Scale(2*block.nodes+block.Flow2Dir(:,end),1));
    block.Effectiveness = FindEffectiveness(Inlet.Flow1,Inlet.Flow2,Flow1,[]);%calculate effectiveness

    %% set up ports : Inlets need to either connected or have initial condition, outlets need an initial condition, and it doesn't matter if they have a connection 
    block.InletPorts = {'Flow1','Flow2','ColdPout','HotPout'};
    block.Flow1.IC = Inlet.Flow1;
    block.Flow1.Saturation = [0,inf];
    block.Flow2.IC = Inlet.Flow2;
    block.Flow2.Saturation = [0,inf];
    block.ColdPout.IC = 101; %Atmospheric pressure
    block.ColdPout.Saturation = [0,inf];
    block.ColdPout.Pstate = []; %identifies the state # of the pressure state if this block has one
    block.HotPout.IC = 101; %Atmospheric pressure
    block.HotPout.Saturation = [0,inf];
    block.HotPout.Pstate = []; %identifies the state # of the pressure state if this block has one
    
    block.OutletPorts = {'ColdOut','HotOut','ColdPin','HotPin'};
    block.ColdOut.IC = Flow1;
    block.HotOut.IC = Flow2;
    block.ColdPin.IC  = block.ColdPout.IC+block.PdropCold;
    block.ColdPin.Pstate = 3*block.nodes+1; %identifies the state # of the pressure state if this block has one
    block.HotPin.IC  = block.HotPout.IC+block.PdropHot;
    block.HotPin.Pstate = 3*block.nodes+2; %identifies the state # of the pressure state if this block has one
    
    block.P_Difference = {'HotPin','HotPout'; 'ColdPin', 'ColdPout';};
    %no dMdP or mFlow (fixed pressure drop)
    Out = block;
elseif length(varargin)==2 %% Have inlets connected, re-initialize
    block = varargin{1};
    Inlet = varargin{2};
    Inlet = checkSaturation(Inlet,block);
    Target = ComponentProperty(block.Target);
    nodes = block.nodes;
    block.ColdPout.IC = Inlet.ColdPout;
    block.HotPout.IC  = Inlet.HotPout;
    block.PfactorCold = NetFlow(Inlet.Flow1)/block.PdropCold;
    block.PfactorHot = NetFlow(Inlet.Flow2)/block.PdropHot;
    block.ColdPin.IC  = Inlet.ColdPout+block.PdropCold;
    block.HotPin.IC  = Inlet.HotPout+block.PdropHot;
    
    specFlow1 = fieldnames(Inlet.Flow1);
    specFlow2 = fieldnames(Inlet.Flow2);
    %% Cold flow
    r = length(block.Flow1Dir(:,1));
    for i = 1:1:length(specFlow1)
        if ~strcmp(specFlow1{i},'T')
            Flow1Out_k.(specFlow1{i})(1:block.nodes,1) = Inlet.Flow1.(specFlow1{i})/r;
        end
    end
    %% Hot flow
    r = length(block.Flow2Dir(:,1));
    for i = 1:1:length(specFlow2)
        if ~strcmp(specFlow2{i},'T')
            Flow2Out_k.(specFlow2{i})(1:block.nodes,1) = Inlet.Flow2.(specFlow2{i})/r;
        end
    end
    AreaOld = block.Area;
    %% adjust effective area to achieve desired effectiveness or temperature
    [block.Area,block.Effectiveness,method] = findArea(Inlet.Flow1,Inlet.Flow2,block.h_conv,Target,block.sizemethod);
    
    %% rescale guessed temperatures to match current inlets and this effectiveness
    Flow1 = Inlet.Flow1;
    Flow1.T  = Inlet.Flow2.T;
    Flow2 = Inlet.Flow2;
    Flow2.T = Inlet.Flow1.T;
    Q1 = enthalpy(Flow1) - enthalpy(Inlet.Flow1);%heat added to the stream
    Q2 = enthalpy(Flow2) - enthalpy(Inlet.Flow2);%heat added to the stream
    if abs(Q1)>abs(Q2)
        Flow2.T = Inlet.Flow2.T - block.Effectiveness*(Inlet.Flow2.T - Flow2.T);
        Q2 = enthalpy(Flow2) - enthalpy(Inlet.Flow2);%heat added to the stream
        ColdEffective = abs(Q2/Q1);
        Flow1.T = Inlet.Flow1.T + ColdEffective*(Flow1.T - Inlet.Flow1.T);
    else
        Flow1.T = Inlet.Flow1.T + block.Effectiveness*(Flow1.T - Inlet.Flow1.T);
        Q1 = enthalpy(Flow1) - enthalpy(Inlet.Flow1);
        HotEffective = abs(Q1/Q2);
        Flow2.T = Inlet.Flow2.T - HotEffective*(Inlet.Flow2.T - Flow1.T);
    end
    Dist1 = (block.Scale(1:nodes) - block.Scale(block.Flow1Dir(1,1)))/(block.Scale(block.Flow1Dir(1,end)) - block.Flow1.IC.T);
    Dist2 = (block.Scale(2*nodes+1:3*nodes) - block.Scale(2*nodes+block.Flow2Dir(1,1)))/(block.Scale(2*nodes+block.Flow2Dir(1,end)) - block.Flow2.IC.T);
    block.Scale(1:nodes) = Inlet.Flow1.T + Dist1*(Flow1.T - Inlet.Flow1.T);
    block.Scale(2*nodes+1:3*nodes) = Inlet.Flow2.T + Dist2*(Flow2.T - Inlet.Flow2.T);
    block.Scale(nodes+1:2*nodes) = (block.Scale(1:nodes) + block.Scale(2*nodes+1:3*nodes))/2;% Average hot and cold side
    %%%
    
    block.HTconv = block.HTconv*block.Area/AreaOld;
    Y = [block.Scale(1:3*nodes);1];
    [T, Y] = ode15s(@(t,y) SolveTempsDynamic(t,y,block,Flow1Out_k,Flow2Out_k,Inlet,method,Target), [0, 1e5], Y);
    Y = Y(end,:)';
    block.Area = Y(end)*block.Area;
    block.HTconv = Y(end)*block.HTconv;

    block.Scale(1:3*block.nodes) = Y(1:end-1);
    
    block.Convection = block.h_conv*block.Area/block.nodes/1000; %all heat transfer coefficients converted to kW/K: thus Q = C*(T1-T2) is in kW
    block.Flow1.IC = Inlet.Flow1;
    block.ColdOut.IC = Inlet.Flow1;
    block.ColdOut.IC.T  = mean(block.Scale(block.Flow1Dir(:,end)));
    block.Flow2.IC = Inlet.Flow2;
    block.HotOut.IC = Inlet.Flow2;
    block.HotOut.IC.T  = mean(block.Scale(2*block.nodes+block.Flow2Dir(:,end)));
    
    [block.Effectiveness,Imbalance] = FindEffectiveness(Inlet.Flow1,Inlet.Flow2,block.ColdOut.IC,block.HotOut.IC);%calculate effectiveness
    Out = block;
else%running the model
    t = varargin{1};
    Y = varargin{2};
    Inlet = varargin{3};
    block = varargin{4};
    string1 = varargin{5};
    
    Inlet = checkSaturation(Inlet,block);
    Flow1 = NetFlow(Inlet.Flow1);
    Flow2 = NetFlow(Inlet.Flow2);

    %% seperate out temperatures
    nodes = block.nodes;
    Tflow1 = Y(1:nodes);
    Tsolid = Y(nodes+1:2*nodes);
    Tflow2 = Y(2*nodes+1:3*nodes);
    Pflow1In = Y(3*nodes+1);
    Pflow2In = Y(3*nodes+2);

    Nflow1Out = block.PfactorCold*(Pflow1In-Inlet.ColdPout);%total cold flow out
    Nflow2Out = block.PfactorHot*(Pflow2In-Inlet.HotPout);%total hot flow out
    specFlow1 = fieldnames(Inlet.Flow1);
    specFlow2 = fieldnames(Inlet.Flow2);

    if strcmp(string1,'Outlet')
        if Nflow1Out<=0
            Nflow1Out = Flow1;
        end
        Flow1Out.T  = mean(Tflow1(block.Flow1Dir(:,end),1));
        for i = 1:1:length(specFlow1)
            if ~strcmp(specFlow1{i},'T')
                Flow1Out.(specFlow1{i}) = Inlet.Flow1.(specFlow1{i})*Nflow1Out/Flow1;
            end
        end
        Flow2Out.T  = mean(Tflow2(block.Flow2Dir(:,end),1));
        for i = 1:1:length(specFlow2)
            if ~strcmp(specFlow2{i},'T')
                Flow2Out.(specFlow2{i}) = Inlet.Flow2.(specFlow2{i})*Nflow2Out/Flow2;
            end
        end
        %% Outlet Ports
        Out.ColdOut = Flow1Out;
        Out.HotOut = Flow2Out;

        Out.ColdPin = Pflow1In;
        Out.HotPin = Pflow2In;
        Tags.(block.name).ColdOut = Tflow1(block.Flow1Dir(:,end),1);
        Tags.(block.name).HotOut = Tflow2(block.Flow2Dir(:,end),1);
        [Tags.(block.name).Effectiveness , Tags.(block.name).NetImbalance] = FindEffectiveness(Inlet.Flow1,Inlet.Flow2,Out.ColdOut,Out.HotOut); %% calculate effectiveness & imbalance
    elseif strcmp(string1,'dY')
        %% Flow1
        Flow1Out_k.T = Tflow1;
        for j = 1:1:length(block.Flow1Dir(1,:));%1:columns
            k = block.Flow1Dir(:,j);
            r = length(k);
            if j==1
                Flow1In_k.T(k,1) = Inlet.Flow1.T;
            else
                Flow1In_k.T(k,1) = Flow1Out_k.T(kprev);
            end
            for i = 1:1:length(specFlow1)
                if ~strcmp(specFlow1{i},'T')
                    Flow1In_k.(specFlow1{i})(k,1) = Inlet.Flow1.(specFlow1{i})/r;
                    Flow1Out_k.(specFlow1{i})(k,1) = Inlet.Flow1.(specFlow1{i})/r;
                end
            end
            kprev = k;
        end

        %% Flow 2
        Flow2Out_k.T = Tflow2;
        for j = 1:1:length(block.Flow2Dir(1,:));%1:columns
            k = block.Flow2Dir(:,j);
            r = length(k);
            if j==1
                Flow2In_k.T(k,1) = Inlet.Flow2.T;
            else
                Flow2In_k.T(k,1) = Flow2Out_k.T(kprev);
            end
            for i = 1:1:length(specFlow2)
                if ~strcmp(specFlow2{i},'T')
                    Flow2In_k.(specFlow2{i})(k,1) = Inlet.Flow2.(specFlow2{i})/r;
                    Flow2Out_k.(specFlow2{i})(k,1) = Inlet.Flow2.(specFlow2{i})/r;
                end
            end
            kprev = k;
        end
        QT = block.HTconv*Y(1:3*nodes) + block.HTcond*Y(1:3*nodes);
        dY = 0*Y;
        %energy flows & sepcific heats
        Hout1 = enthalpy(Flow1Out_k);
        Hin1 = enthalpy(Flow1In_k);
        Cp_1 = SpecHeat(Flow1Out_k);
        Hout2 = enthalpy(Flow2Out_k);
        Hin2 = enthalpy(Flow2In_k);
        Cp_2 = SpecHeat(Flow2Out_k);

        % time constants for states
        tC1 = (block.Vol_Cold*Cp_1(k)*Pflow1In./(block.Ru*Tflow1(k)));
        tC2 = (block.Mass*block.Solid_SpecHeat);
        tC3 = (block.Vol_Hot*Cp_2(k)*Pflow2In./(block.Ru*Tflow2(k)));

        for i=1:1:length(block.Flow1Dir(1,:))
            k = block.Flow1Dir(:,i);
            dY(k)= (QT(k) + Hin1(k) - Hout1(k))./tC1; %Cold flow
            if i>1
                dY(k) = dY(k)+dY(kprev);
            end
            kprev = k;
        end
        dY(nodes+1:2*nodes)= QT(nodes+1:2*nodes)./tC2;  % Solid
        for i=1:1:length(block.Flow2Dir(1,:))
            k = block.Flow2Dir(:,i);
            dY(2*nodes+k)= (QT(2*nodes+k) + Hin2(k) - Hout2(k))./tC3; %Hot flow
            if i>1
                dY(2*nodes+k) = dY(2*nodes+k)+dY(2*nodes+kprev);
            end
            kprev = k;
        end
        n = 3*nodes;
        %% Pressure
        dY(n+1) = (Flow1-Nflow1Out)*block.Ru*Inlet.Flow1.T/(block.Vol_Cold);%working with total flow rates 
        dY(n+2) = (Flow2-Nflow2Out)*block.Ru*Inlet.Flow2.T/(block.Vol_Hot);%working with total flow rates 
        Out = dY;
    end
end
end%Ends function HeatExchanger

function dY = SolveTempsDynamic(t,Y,block,Flow1,Flow2,Inlet,method,Target)
nodes = block.nodes;
block.HTconv = block.HTconv*Y(end);
dY = 0*Y;
QT = block.HTconv*Y(1:3*nodes) + block.HTcond*Y(1:3*nodes);
Flow1.T = Y(1:nodes);
Flow1_In = Flow1;
Flow1_In.T(block.Flow1Dir(:,1),1) = Inlet.Flow1.T;
for j = 1:1:length(block.Flow1Dir(1,:));%1:columns
    k = block.Flow1Dir(:,j);
    if j~=1
        Flow1_In.T(k,1) = Flow1.T(kprev);
    end
    kprev = k;
end

Flow2.T = Y(2*nodes+1:3*nodes);
Flow2_In = Flow2;
Flow2_In.T(block.Flow2Dir(:,1),1) = Inlet.Flow2.T;
for j = 1:1:length(block.Flow2Dir(1,:));%1:columns
    k = block.Flow2Dir(:,j);
    if j~=1
        Flow2_In.T(k,1) = Flow2.T(kprev);
    end
    kprev = k;
end

%energy flows & sepcific heats
Hout1 = enthalpy(Flow1);
Hin1 = enthalpy(Flow1_In);
Hout2 = enthalpy(Flow2);
Hin2 = enthalpy(Flow2_In);
tC = (block.Mass*block.Solid_SpecHeat); %have everything scale with the slower time constant of the plate, until it is initialized

dY(1:nodes)= (QT(1:nodes) + Hin1 - Hout1)./tC; %Cold flow
dY(nodes+1:2*nodes)= QT(nodes+1:2*nodes)/tC;  % Plate
dY(2*nodes+1:3*nodes)= (QT(2*nodes+1:3*nodes) + Hin2 - Hout2)./tC; %Hot flow
if strcmp(method,'ColdT')
    error = (Target - mean(Y(block.Flow1Dir(:,end),1)))/1e5;
elseif strcmp(method,'HotT')
    error = (Target - mean(Y(2*block.nodes+block.Flow2Dir(:,end),1)))/1e5;
elseif strcmp(method,'Effectiveness')
    ColdOut = Inlet.Flow1;
    ColdOut.T =  mean(Y(block.Flow1Dir(:,end),1));
    QT = enthalpy(ColdOut) - enthalpy(Inlet.Flow1);
    Flow1_Max = Inlet.Flow1;
    Flow1_Max.T = Inlet.Flow2.T;
    Flow2_Min = Inlet.Flow2;
    Flow2_Min.T = Inlet.Flow1.T;
    maxQT1 = enthalpy(Flow1_Max) - enthalpy(Inlet.Flow1);
    maxQT2 = enthalpy(Inlet.Flow2) - enthalpy(Flow2_Min);

    error = block.Effectiveness - QT/min(maxQT1,maxQT2);
elseif strcmp(method,'fixed')
    error = 0;
end
dY(end) = error; %slow change in area
end%Ends function solveTempsDynamic

function [Area,Effectiveness,method] = findArea(ColdIn,HotIn,h_conv,Target,method)
%%Ideal heat transfer
QinC = enthalpy(ColdIn);
ColdOut = ColdIn;
ColdOut.T = HotIn.T;
QinH = enthalpy(HotIn);
HotOut = HotIn;
HotOut.T = ColdIn.T;
Qcoldhot = enthalpy(ColdOut) - QinC;
Qhotcold = QinH - enthalpy(HotOut);
if strcmp(method,'ColdT')
    ColdOut.T = Target;
    QoutC = enthalpy(ColdOut);
    Effectiveness = (QoutC - QinC)/min(Qcoldhot,Qhotcold);
    if Effectiveness>.98
        Target=.98;
        method = 'Effectiveness';
    else Qnet = QoutC - QinC;
    end
end
if strcmp(method,'HotT')
    HotOut.T = Target;
    QoutH = enthalpy(HotOut);
    Effectiveness = (QinH - QoutH)/min(Qcoldhot,Qhotcold);
    if Effectiveness>.98
        Target=.98;
        method = 'Effectiveness';
    else Qnet = QinH - QoutH;
    end
end
if strcmp(method,'Effectiveness')
    Effectiveness = Target;
    Qnet = Effectiveness*min(Qcoldhot,Qhotcold);
end
if ~strcmp(method,'ColdT')
    QoutC = QinC+Qnet;
    ColdOut.T = ColdIn.T + Effectiveness*(HotIn.T - ColdOut.T);
    C = NetFlow(ColdOut)*SpecHeat(ColdOut); % Cp*flow
    errorT = 1;
    while abs(errorT)>1e-2
        errorT = (QoutC - enthalpy(ColdOut))/C;
        ColdOut.T = ColdOut.T + errorT;
    end
end
if ~strcmp(method,'HotT')
    QoutH = QinH-Qnet;
    HotOut.T = HotIn.T - Effectiveness*(HotIn.T - ColdOut.T);
    C = NetFlow(HotOut)*SpecHeat(HotOut); % Cp*flow
    errorT = 1;
    while abs(errorT)>1e-2
        errorT = (QoutH - enthalpy(HotOut))/C;
        HotOut.T = HotOut.T + errorT;
    end
end
%find the log mean temperature difference to find Qnet then find the effective area for this deltaT
LMTD = ((HotOut.T - ColdIn.T) - (HotIn.T - ColdOut.T))/(log(HotOut.T - ColdIn.T) - log(HotIn.T - ColdOut.T));
Area = 2*1000*Qnet/(LMTD*h_conv); %effective area in m^2
end %Ends function findArea