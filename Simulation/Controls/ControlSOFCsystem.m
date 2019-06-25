function Out = ControlSOFCsystem(varargin)
% Controls for SOFC system, control blower power, heater bypass, anode recirculation and fuel flow rate
% Tries to control the FC inlet temperature using heater bypass (too hot) and fuel utilization (too cold)
% Tries to control the FC exit temperature using the air flow indirectly through blower power
% Measured voltage determines current and fuel flow to meet power and utilization
% Seven (7) Inlets: Four Targets (temperature set point, temperature differential, steam to carbon ration, and power), T FC exit, T FC inlet, Voltage
% Nine (9) outlets: Four associated with the targets, Heater bypass, blower power, fuel flow rate, current, anode recirculation
% Four (4) states: Target inlet temperature, Heater bypass, target RPM, blower power
% May need to add state for current back in to avoid fuel starvation during step changes
global Tags
F=96485.339; % %Faraday's constant in Coulomb/mole
if length(varargin)==1 %first initialization
    block = varargin{1};
    block.PIdescription = {'Target SOFC Inlet Temperature';'Heater Bypass';'Target Blower Speed';'Blower Power';};
    block.TargetDescription = {'Operating Temperature';'Hot/Cold Temperature Differential';'Steam to Carbon Ratio';'Net Power'};
    Target = zeros(length(block.Target),1);
    for j = 1:1:length(Target)
        Target(j) = ComponentProperty(block.Target{j});
    end
    block.Target  = Target;
    
    block.InletPorts = {};
    block.OutletPorts = {};
    for i = 1:1:length(block.Target)
        block.InletPorts(end+1) = {strcat('Target',num2str(i))};
        block.OutletPorts(end+1) = {strcat('Measured',num2str(i))};
        block.(strcat('Target',num2str(i))).IC = block.Target(i);
        block.(strcat('Target',num2str(i))).Saturation =  [-inf,inf];
        block.(strcat('Measured',num2str(i))).IC = block.Target(i);
    end
    block.InletPorts(end+1:end+4) = {'Hot','Cold','Voltage','BlowerSpeed'};
    block.Hot.IC = ComponentProperty(strcat(block.connections{5},'.IC'));
    block.Hot.Saturation =  [0,inf];
    block.Cold.IC = ComponentProperty(strcat(block.connections{6},'.IC'));
    block.Cold.Saturation =  [0,inf];
    block.Voltage.IC = ComponentProperty(strcat(block.connections{7},'.IC'));
    block.Voltage.Saturation =  [-inf,inf];
    block.BlowerSpeed.IC = ComponentProperty(strcat(block.connections{8},'.IC'));
    block.BlowerSpeed.Saturation =  [0,inf];
    
    block.OutletPorts(end+1:end+5) = {'HeaterBypass','Blower','Current','AnodeRecirc','FuelFlow'};
    block.Cells = ComponentProperty(block.Cells);
    block.Fuel = ComponentProperty(block.Fuel);
    block.Utilization = ComponentProperty(block.Utilization);

    Bypass = ComponentProperty(block.HeaterBypass);
    BlowerPower = ComponentProperty(block.Blower);
    BlowerMass = ComponentProperty('Blower.FlowDesign');
    Speed = block.BlowerSpeed.IC;
    CathodeMass = MassFlow(ComponentProperty('FC1.Flow1Out.IC'));
    ComponentProperty('Blower.FlowDesign',CathodeMass);
    BlowerPower = BlowerPower*CathodeMass/BlowerMass; %scale blower power
    Current = block.Target4.IC*1000/(block.Voltage.IC*block.Cells);
    Recirculation = ComponentProperty(block.AnodeRecirc);
    FuelFlow = (block.Cells*Current/(2*F*block.Utilization*(4*block.Fuel.CH4+block.Fuel.CO+block.Fuel.H2))/1000);

    block.HeaterBypass.IC = Bypass;
    block.Blower.IC = BlowerPower; 
    block.Current.IC =  Current;
    block.AnodeRecirc.IC = Recirculation;
    block.FuelFlow.IC = FuelFlow;
    
    if Recirculation > 0
        %use reciculation from 1st initialization to set 
        Steam2Carbon = block.Target3.IC;
        if Steam2Carbon ==2 %steam to carbon ratio is unaffected by WGS effectiveness
            block.WGSeffective = 0.7;
        else
            block.WGSeffective = effectiveWGS(block.Fuel,FuelFlow,0.7,Steam2Carbon,block.Cells*Current/(2*F*1000),Recirculation);
        end
    end
    
    block.P_Difference = {};
    
    block.IC = [1;Bypass;1;1;]; % inital condition
    block.Scale = [block.Cold.IC;1;Speed;BlowerPower;];
    block.UpperBound = [1.5*block.Cold.IC;1;1.5*Speed;inf];
    block.LowerBound = [0.5*block.Cold.IC;-1;.33*Speed;0];
    Out = block;
elseif length(varargin)==2 %second initialization
    block = varargin{1};
    Inlet = varargin{2};
    Speed = Inlet.BlowerSpeed;
    RPMerror = 0;
    PEN_Temperature = mean(ComponentProperty('FC1.T.Elec'));
    BlowerPower = ComponentProperty('Blower.NominalPower');
    StackPower = Inlet.Target4 + BlowerPower;
    
    block.Current.IC = StackPower*1000/(Inlet.Voltage*block.Cells);
    block.dT_cath_PEN = PEN_Temperature - (mean(Inlet.Hot)+mean(Inlet.Cold))/2; %temperature differencce between cathode and PEN
    TavgError = (Inlet.Target1 - PEN_Temperature)/Inlet.Target2;
    TinletError = 0;
    
    deltaT = (mean(Inlet.Hot)-mean(Inlet.Cold));
    dTerror =(deltaT-(Inlet.Target2))/Inlet.Target2;

    Steam2Carbon = Inlet.Target3;
    
    %% adjust blower mass flow to get deltaT correct
    blowerMassFlow = ComponentProperty('Blower.FlowDesign');
    ComponentProperty('Blower.FlowDesign',blowerMassFlow*(1+.5*dTerror));
    block.Blower.IC = BlowerPower*(1+.5*dTerror);
%     blower = block.Blower.IC
    %%if too cold, add more fuel, it will burn and more heat will be recovered
%     Utilization = block.Utilization*(1-0.04*TavgError)
%     block.Utilization = Utilization;
    
%     %% adjust heat exchanger effectiveness to control inlet temperature
%     HXtarget = ComponentProperty('HX1.Effectiveness');
%     HXtarget = min(0.96,HXtarget*(1 + .05*TavgError));
%     ComponentProperty('HX1.Target',HXtarget);
%     ComponentProperty('HX1.sizemethod','Effectiveness');
    
    %% adjust heat exchanger cold exit temperature to control FC inlet temperature
    HXtarget = ComponentProperty('HX1.Target');
    HXtarget = HXtarget + TavgError*Inlet.Target2;
    ComponentProperty('HX1.Target',HXtarget);

    %%adjust exit temperature to avoid energy feedback during iterations
    Terror = (Inlet.Target1 - PEN_Temperature);
    CathodeOut = ComponentProperty('FC1.Flow2Out.IC');
    CathodeOut.T = CathodeOut.T + Terror;
    ComponentProperty('FC1.Flow1Out.IC',CathodeOut);
    AnodeOut = ComponentProperty('FC1.Flow2Out.IC');
    AnodeOut.T = AnodeOut.T + Terror;
    ComponentProperty('FC1.Flow2Out.IC',AnodeOut);

    FuelFlow = (block.Cells*block.Current.IC/(2*F*block.Utilization*(4*block.Fuel.CH4+block.Fuel.CO+block.Fuel.H2))/1000);
    block.FuelFlow.IC = FuelFlow;
    
    if block.AnodeRecirc.IC > 0
        block.AnodeRecirc.IC = anodeRecircHumidification(block.Fuel,FuelFlow,block.WGSeffective,Steam2Carbon,block.Cells*block.Current.IC/(2*F*1000),block.AnodeRecirc.IC);
    end
    block.Measured1.IC = PEN_Temperature;
    block.Measured2.IC = deltaT;
    block.Measured3.IC = Steam2Carbon;
    block.Measured4.IC = StackPower - BlowerPower;
    
    block.InitializeError = max(abs(TavgError),abs(dTerror)); 
    block.Scale = [mean(Inlet.Cold); 1; Speed; block.Blower.IC;];
    block.IC = [1-TavgError*block.PropGain(1); block.HeaterBypass.IC-TinletError*block.PropGain(2); 1-dTerror*block.PropGain(3); 1-RPMerror*block.PropGain(4);];    
    Out = block;
else
    t = varargin{1};
    Y = varargin{2};
    Inlet = varargin{3};
    block = varargin{4};
    string1 = varargin{5};
    P_Gain = block.PropGain.*block.Scale;
    I_Gain = block.Gain.*block.Scale;
    Speed = Inlet.BlowerSpeed;
    
    averageT = (mean(Inlet.Hot) + mean(Inlet.Cold))/2 + block.dT_cath_PEN;
    TavgError = (Inlet.Target1 - averageT)/Inlet.Target2;
    T_inlet = Y(1)+TavgError*P_Gain(1);
    TinletError = (mean(Inlet.Cold)-T_inlet)/Inlet.Target2; %target a fixed inlet temperature
    Bypass = min(1,max(0,Y(2)+TinletError*P_Gain(2)));
    deltaT = (mean(Inlet.Hot)-T_inlet);%separate the cntrollers by assuming the bypass controller will get the inlet to this target temperature
%     deltaT = (mean(Inlet.Hot)-mean(Inlet.Cold));
    dTerror =(deltaT-Inlet.Target2)/Inlet.Target2; %target a fixed temperature gradient
%     ToutletError = ((mean(Inlet.Hot) - .5*Inlet.Target2 + block.dT_cath_PEN) - Inlet.Target1)/Inlet.Target2; %target a fixed outlet temperature
    TargetSpeed = (Y(3)+dTerror*P_Gain(3));
    SpeedError = (TargetSpeed-Speed)/block.BlowerSpeed.IC ;
    BlowerPower = Y(4)+SpeedError*P_Gain(4);
    StackPower = Inlet.Target4 + BlowerPower;
    Current = StackPower*1000/(Inlet.Voltage*block.Cells);
    Utilization = min(block.Utilization,block.Utilization + .1*(Y(2)+TavgError*P_Gain(2))); %lower utilization if additional pre-heating is necessary
    FuelFlow = block.Cells*Current/(2*F*Utilization*(4*block.Fuel.CH4+block.Fuel.CO+block.Fuel.H2))/1000;

    if block.AnodeRecirc.IC == 0
        Recirculation = 0;
    else
        Steam2Carbon = Inlet.Target3;
        Recirculation = anodeRecircHumidification(block.Fuel,FuelFlow,block.WGSeffective,Steam2Carbon,block.Cells*Current/(2*F*1000),block.AnodeRecirc.IC);
    end

    if strcmp(string1,'Outlet')
        Out.Measured1 = averageT;
        Out.Measured2 = deltaT;
        Out.Measured3 = Steam2Carbon;
        Out.Measured4 = StackPower - BlowerPower;
        Out.HeaterBypass = Bypass;
        Out.Blower = BlowerPower;
        Out.Current = Current;
        Out.AnodeRecirc = Recirculation;
        Out.FuelFlow = FuelFlow;
        Tags.(block.name).Bypass = Bypass;
        Tags.(block.name).Blower = BlowerPower;
        Tags.(block.name).TargetSpeed = TargetSpeed;
        Tags.(block.name).TargetInletT = T_inlet;
        Tags.(block.name).Recirculation = Recirculation;
        Tags.(block.name).FuelFlow = FuelFlow;
        Tags.(block.name).Current = Current;
        Tags.(block.name).Utilization = Utilization;
        Tags.(block.name).Taverage = averageT;
        Tags.(block.name).Tin = mean(Inlet.Cold);
        Tags.(block.name).Tout = mean(Inlet.Hot);
        Tags.(block.name).deltaT = deltaT;
    elseif strcmp(string1,'dY')
        dY = Y*0;
        dY(1) = I_Gain(1)*TavgError;
        dY(2) = I_Gain(2)*TinletError;
        dY(3) = I_Gain(3)*dTerror;
        dY(4) = I_Gain(4)*SpeedError;
        Out = dY;  
    end
end
end%Ends function ControlSOFCsystem