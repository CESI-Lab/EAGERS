function block = InitializeControlFCstack(varargin)
% Controls for Fuel cell stack only, control air flow and inlet temperature and fuel flow rate, net current and anode recirculation
% Four (4) inlets: T hot, T cold, T average, Voltage
% Six (6) outlets: OxidentTemp, oxidant flow rate, fuel temp, fuel flow rate, current, anode recirculation
% Two (2) states: OxidentTemp, oxidant flow rate, 
% if OxyFC: % Three (3) states: anode recirculation, Fuel flow rate, net current
global F
block = varargin{1};
if length(varargin)==1 %first initialization
    if isfield(block,'OxyFC')
        block.description = {'Anode Recirculation';'Fuel Flow';'Fuel Cell Current';};
    else
        block.description = {'Oxidant Inlet Temperature';'Oxidant Flow Rate';'Fuel Cell Current';};
    end
    
    Target = zeros(length(block.Target),1);
    for j = 1:1:length(block.Target)
        Target(j) = ComponentProperty(block.Target{j});
    end
    block.Target = Target;
    
    block.Cells = ComponentProperty(block.Cells);
    block.Fuel = ComponentProperty(block.Fuel);

    Current = sum(ComponentProperty(block.InitConditions{3}));
    if isfield(block,'OxyFC')
        Recirculation = ComponentProperty(block.InitConditions{1});
        FuelFlow = NetFlow(ComponentProperty(block.InitConditions{2}));
        block.Oxidant = ComponentProperty(block.Oxidant);
        OxidantTemp = block.Target(1);
        block.OxidantUtilization = ComponentProperty(block.OxidantUtilization);
        OxidantFlow = Current*block.Cells/(4000*F*block.Oxidant.O2*block.OxidantUtilization);
    else
        Recirculation = ComponentProperty(block.InitialAnodeRecirc);
        block.Utilization = ComponentProperty(block.Utilization);
        FuelFlow = (block.Cells*Current/(2*F*block.Utilization*(4*block.Fuel.CH4+block.Fuel.CO+block.Fuel.H2))/1000);
        OxidantTemp = ComponentProperty(block.InitConditions{1});
        OxidantFlow = NetFlow(ComponentProperty(block.InitConditions{2}));
    end
    
    
    if Recirculation > 0
        %use reciculation from 1st initialization to set 
        Steam2Carbon = block.Target(3);
        if Steam2Carbon ==2 %steam to carbon ratio is unaffected by WGS effectiveness
            block.WGSeffective = 0.7;
        else
            eWGS = .7; %Initial guess of WGS efftiveness
            CH4 = block.Fuel.CH4*FuelFlow;
            r = Recirculation; 
            dWGS = 1e-6;
            error = 1;
            % Inlet = (Inlet + generated - consumed)*r  + New, thus inlet = New/(1-r) + (generated - consumed)*r/(1-r)
            while abs(error)>1e-6
                COin = block.Fuel.CO*FuelFlow/(1-r) + (block.Fuel.CH4 - eWGS*(block.Fuel.CH4+block.Fuel.CO))*FuelFlow*r/(1-r);
                Inlet.FuelMix.H2O = block.Fuel.H2O*FuelFlow/(1-r) + (block.Cells*Current/(2*F*1000) - (block.Fuel.CH4 + (block.Fuel.CH4 + block.Fuel.CO)*eWGS)*FuelFlow)*r/(1-r);
                S2C = Inlet.FuelMix.H2O/(CH4 + 0.5*COin);
                error = Steam2Carbon - S2C;
                eWGS2 = eWGS+dWGS;
                COin2 = block.Fuel.CO*FuelFlow/(1-r) + (block.Fuel.CH4 - eWGS2*(block.Fuel.CH4+block.Fuel.CO))*FuelFlow*r/(1-r);
                H2Oin2 = block.Fuel.H2O*FuelFlow/(1-r) + (block.Cells*Current/(2*F*1000) - (block.Fuel.CH4 + (block.Fuel.CH4 + block.Fuel.CO)*eWGS2)*FuelFlow)*r/(1-r);
                S2C2 = H2Oin2/(CH4 + 0.5*COin2);
                dSdWGS = (S2C2 - S2C)/dWGS;
                if dSdWGS ==0
                    error = 0;
                    eWGS = 1;
                else
                    eWGS = eWGS + error/dSdWGS;
                end
            end
            block.WGSeffective = eWGS;
        end
    end
        
    block.PortNames = {'Hot','Cold','PEN_Temp','Voltage','OxidantTemp','OxidantFlow','AnodeRecirc','FuelFlow','Current'};
    block.Hot.type = 'in';
    block.Hot.IC = block.Target(2)+.5*block.Target(1); 
    block.Cold.type = 'in';
    block.Cold.IC = block.Target(2)-.5*block.Target(1); 
    block.PEN_Temp.type = 'in';
    block.PEN_Temp.IC = block.Target(2);
    block.Voltage.type = 'in';
    block.Voltage.IC = 0.85; %needed for flow outlet
    
    block.OxidantTemp.type = 'out';
    block.OxidantTemp.IC = OxidantTemp;
    block.OxidantFlow.type = 'out';
    block.OxidantFlow.IC = OxidantFlow; 
    block.AnodeRecirc.type = 'out';
    block.AnodeRecirc.IC = Recirculation;
    block.FuelFlow.type = 'out';
    block.FuelFlow.IC = FuelFlow;
    block.Current.type = 'out';
    block.Current.IC = Current; %needed for flow outlet
    

    if isfield(block,'OxyFC')
        block.IC = [Recirculation;1;1;]; % inital condition
        block.Scale = [1;FuelFlow;Current];
    else
        block.IC = [1 1]; % inital condition
        block.Scale =  [OxidantTemp; OxidantFlow;];
    end
    
    block.P_Difference = {};
    
    for i = 1:1:length(block.PortNames)
        if length(block.connections)<i || isempty(block.connections{i})
            block.(block.PortNames{i}).connected={};
        else
            if ischar(block.connections{i})
                block.(block.PortNames{i}).connected = block.connections(i);
            else
                block.(block.PortNames{i}).IC = block.connections{i};
                block.(block.PortNames{i}).connected={};
            end
        end
    end
end

if length(varargin)==2 %% Have inlets connected, re-initialize
    Inlet = varargin{2};

    NetPower = PowerDemandLookup(0);
    Power = block.Current.IC*Inlet.Voltage*block.Cells/1000;
    
    
    averageT = mean(Inlet.PEN_Temp);
    TavgError = (block.Target(1)-averageT)/block.Target(2);
    deltaT = (mean(Inlet.Hot)-mean(Inlet.Cold));
    dTerror =(deltaT-block.Target(2))/block.Target(2);
     
    Steam2Carbon = block.Target(3);
%     dTerror
%     TavgError
    if isfield(block,'OxyFC')
        FuelFlow = block.FuelFlow.IC*(1-.02*TavgError);
        O2flow = block.Cells*block.Current.IC/(4*F*1000)*32*3600*24/1000; %Ton/day
        Parasitic = (1.0101*O2flow^-.202)*(O2flow*1000/24); %parasitic in kW
        PowerError = (NetPower + Parasitic - Power)/NetPower;
        block.Current.IC = block.Current.IC*(1 + PowerError);
%         block.Current.IC = (NetPower + Parasitic)*1000/(Inlet.Voltage*block.Cells);
    else
        block.Current.IC = NetPower*1000/(Inlet.Voltage*block.Cells);
        FuelFlow = (block.Cells*block.Current.IC/(2*F*block.Utilization*(4*block.Fuel.CH4+block.Fuel.CO+block.Fuel.H2))/1000);
        a = .5;
        block.OxidantFlow.IC = block.OxidantFlow.IC*(1+a*dTerror); 
        block.OxidantTemp.IC = block.OxidantTemp.IC + a*(TavgError + .75*dTerror)*block.Target(2);
    end
    
    block.FuelFlow.IC = FuelFlow;
    if block.AnodeRecirc.IC > 0
        r = block.AnodeRecirc.IC; %Initial guess of anode recirculation
        CH4 = block.Fuel.CH4*FuelFlow;
        dr = 1e-6;
        error = 1;
        % Inlet = (Inlet + generated - consumed)*r  + New, thus inlet = New/(1-r) + (generated - consumed)*r/(1-r)
        while abs(error)>1e-6
            COin = block.Fuel.CO*FuelFlow/(1-r) + (block.Fuel.CH4 - block.WGSeffective*(block.Fuel.CH4+block.Fuel.CO))*FuelFlow*r/(1-r);
            Inlet.FuelMix.H2O = block.Fuel.H2O*FuelFlow/(1-r) + (block.Cells*block.Current.IC /(2*F*1000) - (block.Fuel.CH4 + (block.Fuel.CH4 + block.Fuel.CO)*block.WGSeffective)*FuelFlow)*r/(1-r);
            S2C = Inlet.FuelMix.H2O/(CH4 + 0.5*COin);
            error = Steam2Carbon - S2C;
            r2 = r+dr;
            COin2 = block.Fuel.CO*FuelFlow/(1-r2) + (block.Fuel.CH4 - block.WGSeffective*(block.Fuel.CH4+block.Fuel.CO))*FuelFlow*r2/(1-r2);
            H2Oin2 = block.Fuel.H2O*FuelFlow/(1-r2) + (block.Cells*block.Current.IC /(2*F*1000) - (block.Fuel.CH4 + (block.Fuel.CH4 + block.Fuel.CO)*block.WGSeffective)*FuelFlow)*r2/(1-r2);
            S2C2 = H2Oin2/(CH4 + 0.5*COin2);
            dSdr = (S2C2 - S2C)/dr;
            r = r + error/dSdr;
        end
        block.AnodeRecirc.IC = r;
    end
    if isfield(block,'OxyFC')
        block.InitializeError = 0;%max(abs(PowerError));
%         block.Scale = [1, block.FuelFlow.IC];
%         block.IC = [block.AnodeRecirc.IC-dTerror*block.PropGain(1), 1-PowerError*block.PropGain(2)];%-VoltageError*block.PropGain(3)];
        block.Scale = [1, block.FuelFlow.IC, block.Current.IC];
        block.IC = [block.AnodeRecirc.IC-dTerror*block.PropGain(1), 1-PowerError*block.PropGain(2),1];%-VoltageError*block.PropGain(3)];
    else
        block.InitializeError = max(abs(TavgError),abs(dTerror));
        block.Scale = [block.OxidantTemp.IC, block.OxidantFlow.IC];
        block.IC = [1-TavgError*block.PropGain(1),1-dTerror*block.PropGain(2)];
    end
end