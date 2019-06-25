function Out = FuelCell(varargin)
%This function models a fuel cell  it takes input Y states and returns dY: 
%FC/EC model with many states: Temperatures (bi-polar plate, Air/O2, electrolyte, Fuel/Steam, Bi-polar plate if there is an internal reformer), Air side species ( [ CO2, H2O], N2, O2), Fuel/steam side species (CH4, CO, CO2, H2, H2O, N2, O2), [Reformer species] [Rate of internal reforming reactions], Current, flow1 pressure, flow2 pressure
% Five (5) inlets: {'NetCurrent','Flow1','Flow2','Flow2Pout','Flow1Pout'}
% Seven (7) outlets: {'Flow2Out','Flow1Out','Flow2Pin','Flow1Pin','MeasureCurrent','MeasureTpen','MeasureTflow1','MeasureTflow2'}
% Fuel cell current is positive, Flow1 is the cathode cooling, Flow2 is the anode (fuel),  flow 3 is the internal reformer (upstream of the anode)
% Electrolyzer current is negative, Flow1 is the cooling/heating air, Flow2 is the steam input side, flow 3 is the internal methanation (CO2 injection downstream of the anode)
global Tags
if length(varargin)==1 % first initialization
    block = varargin{1};
    block.F = 96485.339; % %Faraday's constant in Coulomb/mole
    block.Ru = 8.314472; % Universal gas constant in kJ/K*kmol
    block.nodes = block.rows*block.columns;

%%%---%%% User defined variables
    %%Electrochemical parameters %%%
    % H2 +1/2O2 --> H2O (Nernst Eo)
    %SOFC Conductivity - Implemented Equation  sigma = A*e^(-deltaG/RT)/T
    switch block.FCtype
        case 'SOFC'
            block.ElecConst = 2e3; %(K/ohm*m) Electrolyte Constant  %default SOFC  = 9e7
            block.deltaG = 8.0e3; %(kJ/kmol)
            %%unused:
%             block.Io = 5000;          % [A/m2] activation current density %default SOFC = 1000
%             block.DeffO2 = 4.0e-5; %(m^2/s)
%             block.alpha=.7;
            block.t_Membrane = 18e-6;                     % [m] thickness of membrane
            block.t_Cath = 800e-6;                        % [m] thickness of cathode structure
            block.t_An = 50e-6;                           % [m] thickness of Anode structure
            block.t_Elec = block.t_Membrane+block.t_Cath+block.t_An;        % [m] thickness of complete electrolyte            

            block.k_Elec =6.19;                                % [W/m K]  Conductivity of the Electrolyte
            block.Density_Elec = 375;                          % [kg/m3] Density of Electrolyte
            block.C_Elec = .800;                                  % [kJ/(kg K)] specific heat of electrolyte 
        case 'MCFC'
            block.Io = 500;            % [A/m2] activation current density
            block.alpha=.4;      
            block.J_L = 6000;            % [A/m2] Limiting  current density   
            block.Cr0 = 4.7833e-4;%4.25e-5;%
            block.Cr1 = -6.6667e-7;%-5e-8;%
            block.t_Elec = 0.003;        % [m] thickness of complete electrolyte
            
            block.k_Elec =6.19;                                % [W/m K]  Conductivity of the Electrolyte
            block.Density_Elec = 375;                          % [kg/m3] Density of Electrolyte
            block.C_Elec = .800;                                  % [kJ/(kg K)] specific heat of electrolyte 
    end
    %%Geometry
    if strcmp(block.Reformer,'internal') || strcmp(block.Reformer,'methanator')
        block = loadGeometry(block,true); %function that loads channel and plate geometry
        block = FlowDir(block,3); %% Load flow direction
    else
        block = loadGeometry(block,false); %function that loads channel and plate geometry
        block = FlowDir(block,2); %% Load flow direction
    end
%%---%%End user defined parameters
    %% Pressure
    block.Flow1_Pout =  block.PressureRatio*101;
    block.Flow1_Pinit = block.Flow1_Pout + block.Flow1Pdrop;
    block.Flow2_Pout = block.PressureRatio*101;
    block.Flow2_Pinit = block.Flow2_Pout + block.Flow2Pdrop;
    
    if strcmp(block.Mode,'fuelcell')
        p_n = 1;
    elseif strcmp(block.Mode,'electrolyzer')
        p_n = -1;
    end
    
    % number of cells
    if strcmp(block.Specification,'cells')
        block.Cells = block.SpecificationValue;
        block.Specification = 'power density';
        block.SpecificationValue = block.RatedStack_kW*100/(block.L_Cell*block.W_Cell*block.Cells);
    elseif strcmp(block.Specification,'power density')
        block.Cells = ceil(abs(block.RatedStack_kW)*100/(block.L_Cell*block.W_Cell*abs(block.SpecificationValue))); %# of cells in stack
    elseif strcmp(block.Specification,'current density')
        block.Cells = ceil(abs(block.RatedStack_kW)*1000/(0.8*1e4*block.L_Cell*block.W_Cell*abs(block.SpecificationValue))); %# of cells in stack (assumes voltage of 0.8)
    elseif strcmp(block.Specification,'voltage')
        block.Cells = ceil(abs(block.RatedStack_kW)*1000/(abs(block.SpecificationValue)*1e4*block.L_Cell*block.W_Cell*0.5)); %# of cells in stack (assumes 0.5 A/cm^2) corrected later
    end 
    %% %% 1st guess at Initial Condition
    
    Current = zeros(block.nodes,1);
    if strcmp(block.Specification,'power density')
        if strcmp(block.Mode,'fuelcell')
            block.Voltage = .85;
        elseif strcmp(block.Mode,'electrolyzer')
            block.Voltage = -1.3;
        end
        i_avg = abs(block.SpecificationValue)/block.Voltage/1000; %convert mW/cm^2 to A/cm^2, assume an initial guess voltage of 0.85
    elseif strcmp(block.Specification,'current density')
        i_avg = p_n*abs(block.SpecificationValue);
        block.Voltage = abs(block.RatedStack_kW)/block.Cells*1000/(block.A_Cell*(100^2))/i_avg; %convert kW to W/cm^2, then divide by A/cm^2 to get V
        Inlet.NetCurrent = i_avg*(block.A_Cell*(100^2));
    elseif strcmp(block.Specification,'voltage')
        block.Voltage = p_n*abs(block.SpecificationValue);
        i_avg = abs(block.RatedStack_kW)/block.Cells*1000/(block.A_Cell*(100^2))/block.Voltage; %convert kW to W/cm^2, then divide by V to get A/cm^2
    end
    for j = 1:1:block.rows
        Current(1+block.columns*(j-1):block.columns*j) =linspace(2,1,block.columns)/sum(linspace(2,1,block.columns))*i_avg*(100^2)*block.A_Cell/block.rows; %make the initial current guess low to not overutilize H2 in 1st iteration of solution
    end
    block.Current.CO = 0*Current;
    block.Current.H2 = Current;
    
    if any(strcmp(block.Reformer,{'external';'adiabatic'})) %% pre-humidify fuel in these cases
        if strcmp(block.Reformer,'external')
            R1 = ComponentProperty('Reformer.ReformTarget')*block.Flow2Spec.CH4;
        elseif strcmp(block.Reformer,'adiabatic')
            R1 = 0.5*block.Flow2Spec.CH4;
        end
        nOut = (1+block.Steam2Carbon*(block.Flow2Spec.CH4+0.5*block.Flow2Spec.CO) + 2*R1);
        ReformedFuel.CH4 = (block.Flow2Spec.CH4 - R1)/nOut;
        ReformedFuel.CO = (block.Flow2Spec.CO + .2*R1)/nOut;
        ReformedFuel.CO2 = (block.Flow2Spec.CO2 + .8*R1)/nOut;
        ReformedFuel.H2 = (block.Flow2Spec.H2 + 3.8*R1)/nOut;
        ReformedFuel.H2O = (block.Flow2Spec.H2O+block.Steam2Carbon*(block.Flow2Spec.CH4+0.5*block.Flow2Spec.CO) - 1.8*R1)/nOut;
        ReformedFuel.N2 = block.Flow2Spec.N2/nOut;
    end
        
    block.T.Flow1 = zeros(block.nodes,1) + block.TpenAvg;
    block.T.Elec = zeros(block.nodes,1) + block.TpenAvg;
    block.T.Flow2 = zeros(block.nodes,1) + block.TpenAvg;
    block.T.Flow3 = zeros(block.nodes,1) + block.TpenAvg;
    %% initial guess of reforming cooling
    if strcmp(block.Mode,'fuelcell')
        FuelSupply  = sum(Current)/(2*block.F*1000)/(block.Utilization_Flow2*(4*block.Flow2Spec.CH4+block.Flow2Spec.CO+block.Flow2Spec.H2))*block.Cells; % fuel flow rate,  current/(2*F*1000) = kmol H2
        if ~isfield(block,'AnPercEquilib')
            block.Recirc.Flow2 = 0;
        else
            block.Recirc.Flow2 = anodeRecircHumidification(block.Flow2Spec,FuelSupply,0.7,block.Steam2Carbon,block.Cells*sum(Current)/(2*block.F*1000),0.5);
        end
        switch block.Reformer
            case 'adiabatic'
                block.ReformT = 823; %an initial guess temperature for adiabatic reforme, or if uncommenting line 873, this is the setpoint
                block.Steam2Carbon = 6; %determines anode recirculation, needs to be high to ensure sufficient temperature for some pre-reforming
        end
        
        switch block.Reformer
            case 'internal'
                block.Flow2_IC  = sum(Current)/(2*block.F*1000)/(block.Utilization_Flow2*(4*block.Flow2Spec.CH4+block.Flow2Spec.CO+block.Flow2Spec.H2))*block.Cells; % fuel flow rate,  current/(2*F*1000) = kmol H2
                R1 = block.Flow2Spec.CH4*block.AnPercEquilib*block.Flow2_IC;
                block.R_CH4ref = block.RefPerc*R1/block.Cells*block.RefSpacing*ones(block.nodes,1)/block.nodes;
                block.R_WGSref =  block.R_CH4ref*.8;
                block.R_CH4 = (R1 - sum(block.R_CH4ref)*block.Cells/block.RefSpacing)/block.Cells*ones(block.nodes,1)/block.nodes;
                block.R_WGS = block.R_CH4*.8;
            case {'external';'adiabatic';}
                block.Flow2Spec = ReformedFuel; %initial fuel composition at inlet
                block.Flow2_IC  = nOut*FuelSupply;%flow out of reformer after recirculation
                R1 = block.Flow2Spec.CH4*block.AnPercEquilib*block.Flow2_IC;
                block.R_CH4 = R1/block.Cells*ones(block.nodes,1)/block.nodes;
                block.R_WGS =  block.R_CH4*.8;
            case {'direct'}
                block.Flow2_IC  = sum(Current)/(2*block.F*1000)/(block.Utilization_Flow2*(4*block.Flow2Spec.CH4+block.Flow2Spec.CO+block.Flow2Spec.H2))*block.Cells; % fuel flow rate,  current/(2*F*1000) = kmol H2
                if isfield(block,'AnPercEquilib')
                    R1 = block.Flow2Spec.CH4*block.AnPercEquilib*block.Flow2_IC;
                else
                    R1 = 0;%no reforming
                end
                block.R_CH4 = R1/block.Cells*ones(block.nodes,1)/block.nodes;
                block.R_WGS =  block.R_CH4*.8;
        end 
        if ~isfield(block,'Utilization_Flow1')
            if block.ClosedCathode
                block.Utilization_Flow1 = 1;
            elseif strcmp(block.Reformer,'internal')
                block.Utilization_Flow1 = .33;
            else
                block.Utilization_Flow1 = .1;
            end
        end
        block.Flow1_IC = block.Cells*sum(Current)/(4*block.F*block.Flow1Spec.O2)/1000/block.Utilization_Flow1;%kmol of oxidant
    elseif strcmp(block.Mode,'electrolyzer')
        block.Flow2_IC  = sum(abs(Current))/(2*block.F*1000)/(block.Utilization_Flow2*block.Flow2Spec.H2O)*block.Cells; % H2O flow rate,  current/(2*block.F*1000) = kmol H2
        block.Recirc.Flow2 = 0;
        block.Flow1_IC = 0;%start with a guess of zero flow on oxidant side
    end
    
    
    
%     %% -- get surface areas and radiation view coefficients from file --%%
%     Dir=strrep(which('InitializeFuelCell.m'),fullfile('Components','Initialization','InitializeFuelCell.m'),'FCMaps');
%     load(fullfile(Dir,block.Map));
%     f = fieldnames(Map);
%     for i = 1:1:length(f)
%         block.(f{i}) = Map.(f{i});
%     end
%     Sigma = 5.670367e-11;%kW/(m^2 K^4) %all heat transfer coefficients converted to kW/m^2*K^4: thus Q = sigma*Area*(T1^4-T2^4) is in kW
%     block.RTmatrix = zeros(s*block.nodes,s*block.nodes);
%     %% Here is where it needs to agregate a 100x100 view factor map into the rows & columns of this particular FC
%     %% -- %%
%     for j = 1:1:block.nodes
%         block.RTmatrix(j,2*block.nodes+1:3*block.nodes) = Sigma*block.ViewFactorCath(j,:)*block.A_Node; %view factor from cathode plate to electrolyte
%         block.RTmatrix(j,j) = -Sigma*sum(block.ViewFactorCath(j,:))*block.A_Node; % - sum(view factors) for this node
%         
%         block.RTmatrix(2*block.nodes+j,1:block.nodes) = Sigma*block.ViewFactorCath(j,:)*block.A_Node; %view factor from electrolyte to cathode plate
%         block.RTmatrix(2*block.nodes+j,2*block.nodes+j) = -Sigma*sum(block.ViewFactorCath(j,:))*block.A_Node; % - sum(view factors) for this node
%         
%         block.RTmatrix(4*block.nodes+j,2*block.nodes+1:3*block.nodes) = Sigma*block.ViewFactorAn(j,:)*block.A_Node; %view factor from anode plate to electrolyte
%         block.RTmatrix(4*block.nodes+j,4*block.nodes+j) = -Sigma*sum(block.ViewFactorAn(j,:))*block.A_Node; % - sum(view factors) for this node
%         
%         block.RTmatrix(2*block.nodes+j,4*block.nodes+1:5*block.nodes) = Sigma*block.ViewFactorAn(j,:)*block.A_Node; %view factor from electrolyte to anode plate
%         block.RTmatrix(2*block.nodes+j,2*block.nodes+j) = block.RTmatrix(2*block.nodes+j,2*block.nodes+j) -Sigma*sum(block.ViewFactorAn(j,:))*block.A_Node; % - sum(view factors) for this node
%         switch block.Reformer
%             case 'internal'
%                 block.RTmatrix(4*block.nodes+j,1:block.nodes) = Sigma*block.ViewFactorRef(j,:)*block.A_Node; %view factor from anode plate to cathode plate, with reformer channels between
%                 block.RTmatrix(4*block.nodes+j,4*block.nodes+j) = block.RTmatrix(4*block.nodes+j,4*block.nodes+j) -Sigma*sum(block.ViewFactorRef(j,:))*block.A_Node; % - sum(view factors) for this node
%                 
%                 block.RTmatrix(j,4*block.nodes+1:5*block.nodes) = Sigma*block.ViewFactorRef(j,:)*block.A_Node; %view factor from cathode plate to anode plate, with reformer channels between
%                 block.RTmatrix(j,j) = block.RTmatrix(j,j) -Sigma*sum(block.ViewFactorRef(j,:))*block.A_Node; % - sum(view factors) for this node
%         end
%     end

   switch block.FCtype
        case 'SOFC'
            criticalSpecies = {'O2';};
        case 'MCFC'
            criticalSpecies = {'CO2';'O2';};
    end
    block.Spec1 = unique([fieldnames(block.Flow1Spec);criticalSpecies]);
    
    criticalSpecies = {'H2';'H2O';};
    if isfield(block.Flow2Spec,'CH4')
        criticalSpecies(end+1) = {'CO';};
        criticalSpecies(end+1) = {'CO2';};
    end
    block.Spec2 = unique([fieldnames(block.Flow2Spec);criticalSpecies]);
    
    Inlet.Flow1.T = block.TpenAvg - .75*block.deltaTStack;
    Inlet.Flow2.T  = block.TpenAvg - block.deltaTStack;
    Inlet.Flow1Pout = block.Flow1_Pout;
    Inlet.Flow2Pout = block.Flow2_Pout;
    
    Inlet = InletFlow(block,Inlet);
    Inlet.Mixed = Inlet.Flow2;%     %%during first initialization it calculates the anode re-cycle and inlet mixing if the fuel is not pre-humidified
    block.T.FuelMix = Inlet.Mixed.T;
    
    %% Run Initial Condition
    [Flow1,Flow2,block,Inlet] = solveInitCond(Inlet,block,1);
    
    if strcmp(block.Reformer,'external') || strcmp(block.Reformer,'adiabatic')% || strcmp(block.Reformer,'none')
        block.Reformer = 'direct'; %external and adiabatic reformers handled in seperate block, after 1st initialization
    end
    %% set up ports : Inlets need to either connected or have initial condition, outlets need an initial condition, and it doesn't matter if they have a connection 
    block.InletPorts = {'NetCurrent','Flow1','Flow2','Flow1Pout','Flow2Pout'};
    block.NetCurrent.IC = sum(block.Current.H2 + block.Current.CO);
    block.NetCurrent.Saturation = [-inf,inf];
    block.Flow1.IC = Inlet.Flow1; 
    block.Flow1.Saturation =  [0,inf];
    block.Flow2.IC = Inlet.Flow2;
    block.Flow2.Saturation =  [0,inf];
    block.Flow1Pout.IC = Inlet.Flow1Pout;
    block.Flow1Pout.Saturation = [0,inf];
    block.Flow1Pout.Pstate = []; %identifies the state # of the pressure state if this block has one
    block.Flow2Pout.IC = Inlet.Flow2Pout;
    block.Flow2Pout.Saturation = [0,inf];
    block.Flow2Pout.Pstate = []; %identifies the state # of the pressure state if this block has one

    block.OutletPorts = {'Flow1Out','Flow2Out','Flow2Pin','Flow1Pin','MeasureVoltage','MeasurePower','MeasureTflow1','MeasureTflow2'};
    block.Flow1Out.IC = MergeLastColumn(Flow1.Outlet,block.Flow1Dir,block.Cells);
    block.Flow2Out.IC  = MergeLastColumn(Flow2.Outlet,block.Flow2Dir,block.Cells);
    block.Flow1Pin.IC = block.Flow1_Pinit;
    block.Flow1Pin.Pstate = length(block.Scale)-1; %identifies the state # of the pressure state if this block has one
    block.Flow2Pin.IC = block.Flow2_Pinit;
    block.Flow2Pin.Pstate = length(block.Scale); %identifies the state # of the pressure state if this block has one
    block.MeasureVoltage.IC = block.Voltage;
    block.MeasurePower.IC = sum((block.Current.H2 + block.Current.CO)*block.Voltage*block.Cells)/1000;%power in kW
    block.MeasureTflow1.IC = block.T.Flow1(block.Flow1Dir(:,end));
    block.MeasureTflow2.IC = block.T.Flow2(block.Flow2Dir(:,end));

    block.P_Difference = {'Flow1Pin', 'Flow1Pout';'Flow2Pin','Flow2Pout';};
    Out = block;
elseif length(varargin)==2 %% Have inlets connected, re-initialize
    block = varargin{1};
    Inlet = varargin{2};
    Inlet = checkSaturation(Inlet,block);
    block.Specification = 'current density';%converge only to match current density from controller
    block.Recirc.Flow2 = 0; %after 1st initialization recirculation is handled in controller, valve & mixing volume
    Inlet.Mixed = Inlet.Flow2; %mixing moved to seperate mixing block
    
    %% Add in species that come from connected blocks that weren't in initial setup
    Flow1New = fieldnames(Inlet.Flow1);
    Flow1All = unique([block.Spec1;Flow1New]);
    Flow1All = Flow1All(~strcmp('T',Flow1All));
    for i = 1:1:length(Flow1All)
        if ~ismember(Flow1All{i},Flow1New)
            Inlet.Flow1.(Flow1All{i})=0;
        end
    end
    block.Spec1 = Flow1All;
    
    Flow2New = fieldnames(Inlet.Flow2);
    Flow2All = unique([block.Spec2;Flow2New]);
    Flow2All = Flow2All(~strcmp('T',Flow2All));
    for i = 1:1:length(Flow2All)
        if ~ismember(Flow2All{i},Flow2New)
            Inlet.Flow2.(Flow2All{i})=0;
        end
    end
    block.Spec2 = Flow2All;
    
    %% Reset pressures based on connected blocks
    block.Flow1_Pinit = Inlet.Flow1Pout + block.Flow1Pdrop;
    block.Flow2_Pinit = Inlet.Flow2Pout + block.Flow2Pdrop;
    block.Flow2Pout.IC = Inlet.Flow2Pout;
    block.Flow1Pout.IC = Inlet.Flow1Pout;
    
    %%--%%
    [Flow1,Flow2,block,~] = solveInitCond(Inlet,block,2);
    %%%
    
    block.Flow1Pin.Pstate = length(block.Scale)-1; %identifies the state # of the pressure state if this block has one
    block.Flow2Pin.Pstate = length(block.Scale); %identifies the state # of the pressure state if this block has one
    block.Flow2Out.IC  = MergeLastColumn(Flow2.Outlet,block.Flow2Dir,block.Cells);
    block.Flow1Out.IC = MergeLastColumn(Flow1.Outlet,block.Flow1Dir,block.Cells);
    block.Flow2Pin.IC = block.Flow2_Pinit;
    block.Flow1Pin.IC = block.Flow1_Pinit;
    block.MeasureCurrent.IC = sum(block.Current.H2 + block.Current.CO);
    block.MeasurePower.IC = abs(sum((block.Current.H2 + block.Current.CO)))*block.Voltage*block.Cells/1000;%power in kW
    block.MeasureTpen.IC = block.T.Elec;
    block.MeasureTflow1.IC = block.T.Flow1(block.Flow1Dir(:,end));
    block.MeasureTflow2.IC = block.T.Flow2(block.Flow2Dir(:,end));
    block.HumidifiedFuelTemp.IC = block.T.FuelMix;
    Out = block;
else%running the model
    t = varargin{1};
    Y = varargin{2};
    Inlet = varargin{3};
    block = varargin{4};
    string1 = varargin{5};
    Inlet = checkSaturation(Inlet,block);
    %% add species that may not be in inlet
    inFields = fieldnames(Inlet.Flow1);
    for i = 1:1:length(block.Spec1)
        if ~ismember(block.Spec1{i},inFields)
            Inlet.Flow1.(block.Spec1{i}) = 0;
        end
    end
    inFields = fieldnames(Inlet.Flow2);
    for i = 1:1:length(block.Spec2)
        if ~ismember(block.Spec2{i},inFields)
            Inlet.Flow2.(block.Spec2{i}) = 0;
        end
    end

    nodes = block.nodes;
    %% seperate out temperatures
    Flow1.Outlet.T = Y(nodes+1:2*nodes);
    T_Elec = Y(2*nodes+1:3*nodes);
    Flow2.Outlet.T = Y(3*nodes+1:4*nodes); 
    switch block.Reformer
        case 'internal'
            nT = 6*nodes; % # of temperature states
            Flow3.Outlet.T = Y(5*nodes+1:6*nodes); %only gets used if internal reformer exists, otherwise these values are actually QT1
        case {'direct';'external';'adiabatic';'none';}
            nT = 4*nodes; % # of temperature states
    end
    n = nT;

    %Current
    nCurrent = Y(end-nodes-1:end-2);%*Inlet.NetCurrent/sum(Y(end-nodes-1:end-2));
    %Pressure
    P_flow1 = Y(end-1); %pressure
    P_flow2 = Y(end); %pressure

    %% Air/Oxidant
    for i = 1:1:length(block.Spec1)
        Flow1.Outlet.(block.Spec1{i}) = Y(n+1:n+nodes); n = n+nodes;
    end
    for j = 1:1:length(block.Flow1Dir(1,:));
        if j==1%first column recieves fresh inlet
            k = block.Flow1Dir(:,1);
            Flow1.Inlet.T(k,1) = Inlet.Flow1.T; 
            for i = 1:1:length(block.Spec1)
                Flow1.Inlet.(block.Spec1{i})(k,1) = Inlet.Flow1.(block.Spec1{i})/block.Cells/length(k); 
            end
        else%subsequent columns recieve outlet of previous column
            k2 = block.Flow1Dir(:,j);
            Flow1.Inlet.T(k2,1) = Flow1.Outlet.T(k);
            for i = 1:1:length(block.Spec1)
                Flow1.Inlet.(block.Spec1{i})(k2,1) = Flow1.Outlet.(block.Spec1{i})(k); 
            end
            k = k2;
        end
        if block.ClosedCathode %closed end cathode
            Flow1.Outlet.O2(k,1) = Flow1.Inlet.O2(k,1) - nCurrent(k)/(4*block.F*1000); 
            if any(strcmp(block.FCtype,{'MCFC','MCEC'}))
                Flow1.Outlet.CO2(k,1) = Flow1.Inlet.CO2(k,1) - nCurrent(k)/(2*block.F*1000); 
            end
        end
    end

    %% Fuel/Steam
    for i = 1:1:length(block.Spec2)
        Flow2.Outlet.(block.Spec2{i}) = Y(n+1:n+nodes); n = n+nodes;
    end
    %% Internal reformer
    switch block.Reformer
        case 'internal'
            for i = 1:1:length(block.Spec2)
                Flow3.Outlet.(block.Spec2{i}) = Y(n+1:n+nodes); n = n+nodes;
            end
    end

    Flow1Out.T  = mean(Flow1.Outlet.T(block.Flow1Dir(:,end))); %temperature 
    for i = 1:1:length(block.Spec1)
        Flow1Out.(block.Spec1{i}) = max(0,sum(Flow1.Outlet.(block.Spec1{i})(block.Flow1Dir(:,end)))*block.Cells);%avoid sending negative outlets
    end
    Flow2Out.T  = mean(Flow2.Outlet.T(block.Flow2Dir(:,end))); %temperature 
    for i = 1:1:length(block.Spec2)
        Flow2Out.(block.Spec2{i}) = max(0,sum(Flow2.Outlet.(block.Spec2{i})(block.Flow2Dir(:,end)))*block.Cells);%avoid sending negative outlets
    end
    %% Reformer
    switch block.Reformer
        case 'internal'
            for j = 1:1:length(block.Flow3Dir(1,:))
                if j==1
                    k = block.Flow3Dir(:,1);
                    Flow3.Inlet.T(k,1) = Inlet.Flow2.T;
                    for i = 1:1:length(block.Spec2)
                        Flow3.Inlet.(block.Spec2{i})(k,1) = Inlet.Flow2.(block.Spec2{i})/block.Cells/length(k)*block.RefSpacing;
                    end
                else
                    k2 = block.Flow3Dir(:,j);
                    Flow3.Inlet.T(k2,1) = Flow3.Outlet.T(k);
                    for i = 1:1:length(block.Spec2)
                        Flow3.Inlet.(block.Spec2{i})(k2,1) = Flow3.Outlet.(block.Spec2{i})(k);
                    end
                    k = k2;
                end
            end
    end
    %% Fuel/Steam
    for j = 1:1:length(block.Flow2Dir(1,:))
        k2 =block.Flow2Dir(:,j);
        if j==1 % first column of fuel flow direction
            switch block.Reformer
                case 'internal'
                    Flow2.Inlet.T(k2,1) =  Flow3.Outlet.T(k);
                case {'direct';'external';'adiabatic';'none';}
                    Flow2.Inlet.T(k2,1) = Inlet.Flow2.T;
            end
            for i = 1:1:length(block.Spec2)
                switch block.Reformer
                    case 'internal'
                        Flow2.Inlet.(block.Spec2{i})(k2,1) = Flow3.Outlet.(block.Spec2{i})(k)/block.RefSpacing;%Species flows coming into anode
                    case {'direct';'external';'adiabatic';'none';}
                        Flow2.Inlet.(block.Spec2{i})(k2,1) = Inlet.Flow2.(block.Spec2{i})/block.Cells/length(k2);
                end
            end
        else
            Flow2.Inlet.T(k2,1) = Flow2.Outlet.T(k);
            for i = 1:1:length(block.Spec2)
                Flow2.Inlet.(block.Spec2{i})(k2,1) = Flow2.Outlet.(block.Spec2{i})(k);
            end
        end
        k = k2;
    end

    %% Nernst & Losses
    FuelCellNernst(Flow1,Flow2,nCurrent,T_Elec,P_flow1,block)
    Voltage = sum(Tags.(block.name).nVoltage'.*(nCurrent/sum(nCurrent)));
    Current.CO = Tags.(block.name).I_CO';
    Current.H2 = Tags.(block.name).I_H2';
    [h,hs] = enthalpy(T_Elec,{'H2','H2O','O2','CO','CO2'});
    nPower = abs(Voltage)*nCurrent/1000; %cell power in kW
    Qreaction = Current.H2/(2000*block.F).*(h.H2+.5*h.O2-h.H2O) + Current.CO/(2000*block.F).*(h.CO+.5*h.O2-h.CO2);
    Qgen = Qreaction - nPower;%kW of heat generated by electrochemistry (per node & per cell) (will be positive in both electrolyzer and fuel cell modes)
    Power = Voltage*abs(Inlet.NetCurrent)*block.Cells/1000;
    if strcmp(string1,'Outlet')
        %% Outlet Ports
        Out.Flow1Out = Flow1Out;
        Out.Flow2Out  = Flow2Out;
        Out.Flow1Pin = P_flow1;
        Out.Flow2Pin = P_flow2;
        Out.MeasureVoltage = Voltage;
        Out.MeasurePower = Power;
        Out.MeasureTflow1 = mean(Y(nodes+block.Flow1Dir(:,end)));
        Out.MeasureTflow2 = mean(Y(3*nodes+block.Flow2Dir(:,end)));
        %% Tags
        Tags.(block.name).Voltage = Voltage;
        Tags.(block.name).Q_gen = sum(Qgen*block.Cells); %kW of heat generated by electrochemistry
        if Voltage>0
            H2_in = sum(Flow2.Inlet.H2(block.Flow2Dir(:,1)));
            H2_out = sum(Flow2.Outlet.H2(block.Flow2Dir(:,end)));
            if isfield(Flow2.Inlet,'CH4')
                H2_in = H2_in + 4*sum(Flow2.Inlet.CH4(block.Flow2Dir(:,1)));
                H2_out = H2_out + 4*sum(Flow2.Outlet.CH4(block.Flow2Dir(:,end)));
            end
            if isfield(Flow2.Inlet,'CO')
                H2_in = H2_in + sum(Flow2.Inlet.CO(block.Flow2Dir(:,1)));
                H2_out = H2_out + sum(Flow2.Outlet.CO(block.Flow2Dir(:,end)));
            end
            Tags.(block.name).H2utilization = (H2_in - H2_out)./H2_in;
            Tags.(block.name).H2concentration = (Flow2.Inlet.H2 + Flow2.Outlet.H2)./(NetFlow(Flow2.Inlet)+NetFlow(Flow2.Outlet));
            Tags.(block.name).O2utilization = sum(Flow1.Inlet.O2(block.Flow1Dir(:,1)) - Flow1.Outlet.O2(block.Flow1Dir(:,end)))/sum(Flow1.Inlet.O2(block.Flow1Dir(:,1)));
            if strcmp(block.FCtype,'MCFC')
                Tags.(block.name).CO2utilization = sum(Flow1.Inlet.CO2(block.Flow1Dir(:,1)) - Flow1.Outlet.CO2(block.Flow1Dir(:,end)))/sum(Flow1.Inlet.CO2(block.Flow1Dir(:,1))); %only makes sense if FCtype=1 and CO2 is a cathode state
            else
                Tags.(block.name).CO2utilization =0;
            end
            Tags.(block.name).H2Outilization = 0;
        else
            Tags.(block.name).H2Outilization = sum(Flow2.Inlet.H2O(block.Flow2Dir(:,1)) - Flow2.Outlet.H2O(block.Flow2Dir(:,end)))/sum(Flow2.Inlet.H2O(block.Flow2Dir(:,1)));
            Tags.(block.name).H2utilization = 0;
            Tags.(block.name).O2utilization = 0;
            Tags.(block.name).CO2utilization = 0;
        end
        Tags.(block.name).T_flow1_out = Flow1Out.T;
        Tags.(block.name).T_flow2_out = Flow2Out.T;
        Tags.(block.name).Current = sum(nCurrent);
        Tags.(block.name).StackPower = Tags.(block.name).Current*abs(Voltage)*block.Cells/1000; %power in kW
        Tags.(block.name).StackdeltaT = Flow1Out.T-Inlet.Flow1.T;
        if isfield(block,'CoolingStream') && strcmp(block.CoolingStream,'anode')
            Tags.(block.name).StackdeltaT = Flow2Out.T-mean(Flow2.Inlet.T(block.Flow2Dir(:,1)));
        end
        Tags.(block.name).PENtemperature = T_Elec';
        Tags.(block.name).PENavgT = sum(T_Elec)/block.nodes;
        Tags.(block.name).MaxPEN = max(T_Elec);
        Tags.(block.name).PENdeltaT = Tags.(block.name).MaxPEN-min(T_Elec);
        Tags.(block.name).dTdX = (T_Elec-T_Elec(block.HTadjacent(:,2)))/(block.L_Cell/block.columns);
        Tags.(block.name).dTdY = (T_Elec-T_Elec(block.HTadjacent(:,4)))/(block.W_Cell/block.rows);
        Tags.(block.name).MaxdTdX = max(abs([Tags.(block.name).dTdX;Tags.(block.name).dTdY;]));
        Tags.(block.name).Efficiency = Power/(NetFlow(Inlet.Flow2)*HeatingValue(Inlet.Flow2));
        
    elseif strcmp(string1,'dY')
        if ~isfield(block,'AnPercEquilib') %no reforming or WGS to worry about
            spec = fieldnames(Flow2.Inlet);
            spec = spec(~strcmp('T',spec));
            for i = 1:1:length(spec)
                if strcmp(spec{i},'CO2') 
                    Flow2Out.CO2 = Flow2.Inlet.CO2 + Current.CO/(2*block.F*1000); 
                    if strcmp(block.FCtype,'MCFC')
                        Flow2Out.CO2 = Flow2Out.CO2 + (Current.H2+Current.CO)/(2*block.F*1000);
                    end
                elseif strcmp(spec{i},'H2') 
                    Flow2Out.H2 = (Flow2.Inlet.H2 - Current.H2/(2*block.F*1000));
                elseif strcmp(spec{i},'H2O') 
                    Flow2Out.H2O = (Flow2.Inlet.H2O + Current.H2/(2*block.F*1000)); 
                elseif strcmp(spec{i},'CO') 
                    Flow2Out.CO = (Flow2.Inlet.CO - Current.CO/(2*block.F*1000));  
                else
                    Flow2Out.(spec{i}) = Flow2.Inlet.(spec{i});  
                end
            end
        else
            [R,Flow2Out] = KineticReformation(block.method,Flow2,P_flow2,block.KineticCoeff1,Current,block);%% Kinetic reaction rates  (WGS is always near equilibrium)
        end
        switch block.Reformer
            case 'internal'
                RefCurrent.H2 = zeros(nodes,1);
                RefCurrent.CO = zeros(nodes,1);
                [Rref,RefOut] = KineticReformation(block.method,Flow3,P_flow2,block.KineticCoeff3,RefCurrent,block);%% Kinetic reaction rates  (WGS is always near equilibrium)
                R.CH4ref = Rref.CH4;
                R.WGSref = Rref.WGS;
        end
        
        switch block.FCtype%ion transport across membrane (total enthalpy), 
            case {'SOFC';}%ion crosses from flow 1 to flow 2 in fuel cell mode and from flow 2 to flow 1 in electrolyzer mode. Current is negative in electrolyzer mode
                Qion = (nCurrent)/(4000*block.F).*hs.O2; %O2 ion crossing over (kW)
            case {'MCFC';}%ion crosses from flow 1 to flow 2 in fuel cell mode and from flow 2 to flow 1 in electrolyzer mode. Current is negative in electrolyzer mode
                Qion = (nCurrent)/(4000*block.F).*hs.O2 + abs(nCurrent)/(2000*block.F).*hs.CO2;% O2 & CO2 ion crossing over
        end
        %% Q %% Heat transfer & Generation
        switch block.Reformer
            case 'internal'
                QT = block.HTconv*Y(1:6*nodes) + block.HTcond*Y(1:6*nodes) + block.HTrad*(Y(1:6*nodes).^4);
            case {'direct';'external';'adiabatic';'none';}
                QT = block.HTconv*Y(1:4*nodes) + block.HTcond*Y(1:4*nodes) + block.HTrad*(Y(1:4*nodes).^4);
        end

        %energy flows & sepcific heats
        Hout1 = enthalpy(Flow1.Outlet);
        Hin1 = enthalpy(Flow1.Inlet);
        Hout2 = enthalpy(Flow2.Outlet);
        Hin2 = enthalpy(Flow2.Inlet);

        %% %% solve for dY in order of states
        dY = 0*Y;
        %%Temperatures
        dY(1:nodes)= QT(1:nodes)./block.tC(1:nodes);  %Bi-Polar Plate
        for i=1:1:length(block.Flow1Dir(1,:)) %having the downstream nodes change temperature with the upstream nodes prevents propogation issues when taking larger time steps
            k = block.Flow1Dir(:,i);
            dY(nodes+k)= (QT(nodes+k) + Hin1(k) - Hout1(k) - Qion(k))./block.tC(nodes+k); %Air/O2: ion always leaves this node
            if i>1
                dY(nodes+k) = dY(nodes+k)+dY(nodes+kprev);
            end
            kprev = k;
        end
        dY(1+2*nodes:3*nodes)= (QT(1+2*nodes:3*nodes) + Qgen)./block.tC(2*nodes+1:3*nodes); %Electrolyte Plate
        for i=1:1:length(block.Flow2Dir(1,:))
            k = block.Flow2Dir(:,i);
            dY(3*nodes+k)= (QT(3*nodes+k) + Hin2(k) - Hout2(k) + Qion(k) - Qreaction(k))./block.tC(3*nodes+k);  %fuel/steam
            if i>1
                dY(3*nodes+k) = dY(3*nodes+k)+dY(3*nodes+kprev);
            end
            kprev = k;
        end
        n =nT;
        %%Air/oxidant Species
        for i = 1:1:length(block.Spec1)
            if strcmp(block.Spec1{i},'CO2') && any(strcmp(block.FCtype,{'MCFC';'MCEC'}))
                dY(n+1:n+nodes)= (Flow1.Inlet.CO2 - Y(n+1:n+nodes) - nCurrent/(2*block.F*1000))./block.tC(n+1:n+nodes);  %CO2 species concentration with CO2 crossover
            elseif strcmp(block.Spec1{i},'O2') && any(strcmp(block.FCtype,{'SOFC';'SOEC';'MCFC';'MCEC'}))
                dY(n+1:n+nodes)= (Flow1.Inlet.O2 - Y(n+1:n+nodes) - nCurrent/(4*block.F*1000))./block.tC(n+1:n+nodes);%O2 species concentration with O2 crossover
            else
                dY(n+1:n+nodes)= (Flow1.Inlet.(block.Spec1{i}) - Flow1.Outlet.(block.Spec1{i}))./block.tC(n+1:n+nodes);%all other species concentration
            end
            n = n+nodes;
        end
        %% Fuel/steam Species
        for i = 1:1:length(block.Spec2)
            dY(n+1:n+nodes)= (Flow2Out.(block.Spec2{i}) - Flow2.Outlet.(block.Spec2{i}))./block.tC(n+1:n+nodes); %all species concentration
            n = n+nodes; 
        end
        %% Reformer
        switch block.Reformer
            case 'internal'
                dY(1+4*nodes:5*nodes)= QT(1+4*nodes:5*nodes)./block.tC(4*nodes+1:5*nodes);  %2nd half of Plate
                
                Hout3 = enthalpy(Flow3.Outlet);
                Hin3 = enthalpy(Flow3.Inlet);
                for i=1:1:length(block.Flow3Dir(1,:))
                    k = block.Flow3Dir(:,i);
                    dY(5*nodes+k)= (block.RefSpacing*QT(5*nodes+k) + Hin3(k) - Hout3(k))./block.tC(5*nodes+k);  %Fuel Reformer Channels
                    if i>1
                        dY(5*nodes+k) = dY(5*nodes+k)+dY(5*nodes+kprev);
                    end
                    kprev = k;
                end
                for i = 1:1:length(block.Spec2)
                    dY(n+1:n+nodes)= (RefOut.(block.Spec2{i}) - Flow3.Outlet.(block.Spec2{i}))./block.tC(n+1:n+nodes);   %all species concentrations
                    n = n+nodes;
                end
        end

        %%Current
        dY(n+1:n+nodes) = (Inlet.NetCurrent - sum(Y(end-nodes-1:end-2)) + (abs(Tags.(block.name).nVoltage)'-abs(Voltage))./Tags.(block.name).ASR.*(block.A_Cell*100^2))./block.tC(end-nodes-1:end-2); n = n+nodes; %error in A/cm^2 * area  

        %%Pressure
        if block.ClosedCathode %closed end cathode
            dY(n+1) = 0; %no flow out, so anything not used in electrochemistry adds to pressure
        else
            Nflow1 = block.Pfactor1*max(0.01,(P_flow1-Inlet.Flow1Pout));%total air/oxidant flow out
            dY(n+1) = (NetFlow(Inlet.Flow1)-Nflow1)*block.Ru*Inlet.Flow1.T/block.tC(n+1);
        end
        Nflow2 = block.Pfactor2*max(0.1,(P_flow2-Inlet.Flow2Pout));%total fuel/steam flow out
        dY(n+2) = (NetFlow(Inlet.Flow2)-Nflow2)*block.Ru*Inlet.Flow2.T/block.tC(n+2);%working with total flow rates so must multiply by nodes & cells
        Out = dY;
    end
end
end %Ends function FuelCell

function [Flow1,Flow2,block,Inlet] = solveInitCond(Inlet,block,firstSolve)
error = 1;
Tol = 1e-3;
count = 1;
while abs(error)>Tol %iterate to reach target current density, voltage or power density
    Flow1 = FCin2Out(block.T.Flow1,Inlet.Flow1,block.Flow1Dir,block.FCtype,block.Cells,block.Current,[],'Flow1');
%     SinglePassUtilization = (sum(block.Current)*block.Cells/(2000*F))/(4*Inlet.Mixed.CH4+Inlet.Mixed.CO + Inlet.Mixed.H2);
    if ~isfield(block,'AnPercEquilib') %no reforming or WGS to worry about
        Flow2 = FCin2Out(block.T.Flow2,Inlet.Flow2,block.Flow2Dir,block.FCtype,block.Cells,block.Current,[],'Flow2');
        Flow3 = [];%no reformer
    else
        [block,Inlet.Mixed,Flow2,Flow3] = KineticCoef(block,Inlet,(firstSolve==1),(count==1));%% solve for kinetic reaction coefficient which results in the correct AnPercEquilib & RefPerc (match R_CH4 & R_WGS)
    end
    Offset = 0;
    if count==1 && firstSolve==1
        F1spec = block.Flow1Spec;%neccessary when the flow rate is 0
        F1spec.T = Inlet.Flow1.T;
        O2consume = block.Cells*sum(block.Current.H2 + block.Current.CO)/(4000*block.F);
        mdot_Cp(1) = SpecHeat(F1spec).*(NetFlow(Inlet.Flow1)-0.5*O2consume)/(block.Cells*block.rows);
        
        if strcmp(block.direction ,'crossflow')
            mdot_Cp(2) = SpecHeat(Inlet.Flow2).*NetFlow(Inlet.Flow2)/(block.Cells*block.columns);
        else
            mdot_Cp(2) = SpecHeat(Inlet.Flow2).*NetFlow(Inlet.Flow2)/(block.Cells*block.rows);
        end
        [block.Tstates,block.HTcond,block.HTconv,block.HTrad]= SteadyTemps(block,mdot_Cp,[Inlet.Flow1.T,Inlet.Mixed.T]);
        block = Set_IC(block,Flow1,Flow2,Flow3);%do this to get time constants tC
    else
        [~, Y] = ode15s(@(t,y) DynamicTemps(t,y,block,Flow1,Flow2,Flow3,Inlet), [0, 1e4], block.Tstates);
        block.Tstates = Y(end,:)';
        if firstSolve==1 && abs(mean(block.Tstates(2*block.nodes+1:3*block.nodes))-block.TpenAvg)>10 %temperature in solve dynamic is diverging too much (1300K), and messing up reforming solution (this is a temporary fix
            Offset = (mean(block.Tstates(2*block.nodes+1:3*block.nodes))-block.TpenAvg);
        end
    end
    %organize temperatures
    block.T.Flow1 =  block.Tstates(1*block.nodes+1:2*block.nodes) - Offset;
    block.T.Elec =  block.Tstates(2*block.nodes+1:3*block.nodes) - Offset;
    block.T.Flow2 =  block.Tstates(3*block.nodes+1:4*block.nodes) - Offset;
    block.T.FuelMix = Inlet.Mixed.T;
    switch block.Reformer
        case 'internal'
            block.T.Flow3 = block.Tstates(5*block.nodes+1:6*block.nodes) - Offset;
    end
    T = block.TpenAvg + (block.T.Elec-mean(block.T.Elec)); %assume you will get to the desired temperature (this avoids oscilations in voltage and helps convergence
    FuelCellNernst(Flow1,Flow2,block.Current,T,block.Flow1_Pinit,block);
    Current(count) = {abs(sum(block.Current.H2 + block.Current.CO))};
    Voltage(count) = {block.Voltage};
    [block,error,scale] = redistributeCurrent(block,Inlet,Current{max(1,count-1)},Voltage{max(1,count-1)},count,firstSolve);%% calculate the change in current to converge to the desired power density, voltage, or current
    TotCurrent = abs(sum(block.Current.H2 + block.Current.CO));
    
    if strcmp(block.Mode,'fuelcell') && block.ClosedCathode % Ensure there is enough oxidant so Flow2 does not go negative
        Inlet.Flow1.O2 = block.Cells*TotCurrent/(4000*block.F*block.Flow1Spec.O2)/block.Utilization_Flow1;%kmol of oxidant
        if any(strcmp(block.FCtype,{'MCFC'}))
            Inlet.Flow1.CO2 = block.Cells*TotCurrent/(2000*block.F*block.Flow1Spec.O2)/block.Utilization_Flow1;%kmol of oxidant
        end
    end
    if firstSolve ==1 %solving block to convergence without other blocks or controller
        if isfield(block,'R_CH4')
            block.R_CH4 = scale*block.R_CH4;
            block.R_WGS = scale*block.R_WGS;
        end
        if isfield(block,'R_CH4ref')
            block.R_CH4ref = scale*block.R_CH4ref;% fuel flow scales with current so assume reforming will
            block.R_WGSref = scale*block.R_WGSref;
        end
        if strcmp(block.Mode,'fuelcell')
            if strcmp(block.CoolingStream,'cathode')%air flow balances heat generation to maintain deltaT, heat transfer to anode and any fuel reforming is accounted for
                block.Flow2_IC  = TotCurrent/(2*block.F*1000)/(block.Utilization_Flow2*(4*block.Flow2Spec.CH4+block.Flow2Spec.CO+block.Flow2Spec.H2))*block.Cells; % Fresh fuel flow rate,  current/(2*F*1000) = kmol H2
                k = block.Flow1Dir(:,end);
                dTerror = (mean((block.Tstates(k+block.nodes))- Inlet.Flow1.T)/block.deltaTStack-1);
                block.Flow1_IC = block.Flow1_IC*(1 + dTerror)*scale^2;
                TavgError = (block.TpenAvg-mean(block.Tstates(2*block.nodes+1:3*block.nodes)))/block.deltaTStack;
                Inlet.Flow1.T = Inlet.Flow1.T + (TavgError + 0.75*dTerror)*block.deltaTStack;
            elseif strcmp(block.CoolingStream,'anode') %oxidant flow rate determined by current, fuel flow rate is now determined by thermal balancing
                block.Flow1_IC = block.Cells*TotCurrent/(4000*block.F*block.Flow1Spec.O2)/block.Utilization_Flow1;%kmol of oxidant
                if block.ClosedCathode %%energy balance
                    Hin1 = enthalpy(Flow1.Inlet);
                    Hout1 = enthalpy(Flow1.Outlet);
                    Hin2 = enthalpy(Flow2.Inlet);
                    Hout2 = enthalpy(Flow2.Outlet);
                    Hin3 = enthalpy(Flow3.Inlet);
                    Hout3 = enthalpy(Flow3.Outlet);
                    Power = abs(block.Voltage)*(block.Current.H2 + block.Current.CO)/1000; %cell power in kW
                    Qimbalance = sum((Hin2 - Hout2) + (Hin1 - Hout1) + (Hin3 - Hout3) - Power);
                    h = enthalpy(mean(block.T.Flow3),{'H2','H2O','O2','CO','CO2','CH4'});
                    Qreform = (h.CO+3*h.H2-h.CH4-h.H2O) + 0.8*(h.CO2+h.H2-h.CO-h.H2O); %kW of cooling per kmol of fuel
                    ExtraFuel = 0.75*Qimbalance*block.Cells/Qreform/block.Flow2Spec.CH4;
                    error = max(abs(error),abs(Qimbalance/sum(Power)));
                else
                    ExtraFuel = 0;
                    %need to do something with recirculation
                end
                block.Flow2_IC  = block.Cells*TotCurrent/(2000*block.F)/(block.Utilization_Flow2*(4*block.Flow2Spec.CH4+block.Flow2Spec.CO+block.Flow2Spec.H2)); % re-calculate with revised current
                block.Flow2_IC = block.Flow2_IC + ExtraFuel;
                block.Utilization_Flow2 = block.Cells*TotCurrent/(2000*block.F)/(block.Flow2_IC*(4*block.Flow2Spec.CH4+block.Flow2Spec.CO+block.Flow2Spec.H2));
                % change steam to carbon to affect deltaT?
            end
        elseif strcmp(block.Mode,'electrolyzer')
            block.Flow2_IC  = block.Cells*TotCurrent/(2000*block.F*block.Utilization_Flow2*block.Flow2Spec.H2O); % Fresh steam flow rate,  current/(2*F*1000) = kmol H2O
            Hout2 = enthalpy(Flow2.Outlet);
            Hin2 = enthalpy(Flow2.Inlet);
            Q_steam = block.Cells*sum(Hout2(block.Flow2Dir(:,end)) - Hin2(block.Flow2Dir(:,1)));
            Q2oxidant = block.Cells*block.Voltage*TotCurrent/1000 + Q_steam;
            block.Flow1_IC = abs(Q2oxidant)/(33*block.deltaTStack); %air flow is extra heat / Cp* deltaT
            Inlet.Flow1.T = Inlet.Flow1.T + (block.TpenAvg-mean(block.Tstates(2*block.nodes+1:3*block.nodes)))/block.deltaTStack;
        end
        Inlet = InletFlow(block,Inlet);
    end
    count= count+1;
end
Flow1 = FCin2Out(block.T.Flow1,Inlet.Flow1,block.Flow1Dir, block.FCtype,block.Cells,block.Current,[],'Flow1');
block.Pfactor1 = NetFlow(Inlet.Flow1)/block.Flow1Pdrop;
block.Pfactor2 = NetFlow(Inlet.Flow2)/block.Flow2Pdrop;
block = Set_IC(block,Flow1,Flow2,Flow3);
end %Ends function solveInitCond

function Inlet = InletFlow(block,Inlet) %only used 1st time through initialization (before we know what is connected to inlet
%Air/oxidant
switch block.FCtype
    case 'SOFC'
        if block.ClosedCathode
            Inlet.Flow1.O2 = block.Flow1Spec.O2*block.Flow1_IC;
        else
            for i = 1:1:length(block.Spec1)
                if isfield(block.Flow1Spec,block.Spec1{i})
                    Inlet.Flow1.(block.Spec1{i}) = block.Flow1Spec.(block.Spec1{i})*block.Flow1_IC;
                else Inlet.Flow1.(block.Spec1{i}) = 0;
                end
            end
        end
    case 'MCFC' %recalculate cathode inlet species for MCFC (this is an estimate assuming the 100% of non-recirculated anode gas is oxidized and fed to the cathode)
        Inlet.Flow2.CO2 = (Inlet.Mixed.CH4+Inlet.Mixed.CO+Inlet.Mixed.CO2) + sum(block.Current.H2+block.Current.CO)/(2*F*1000)*block.Cells;
        Inlet.Flow2.H2O = (4*Inlet.Mixed.CH4+Inlet.Mixed.CO+Inlet.Mixed.H2+Inlet.Mixed.H2O);
        nonCO2_H2O = (block.Flow2_IC - Inlet.Flow2.CO2 - Inlet.Flow2.H2O);
        for i = 1:1:length(block.Spec2)
            if isfield(block.Flow2Spec,block.Spec2{i})
                if strcmp(block.Spec2{i},'CO2')||strcmp(block.Spec2{i},'H2O')
                    Inlet.Flow2.(block.Spec2{i}) = Inlet.Flow2.(block.Spec2{i}) + block.Flow2Spec.(block.Spec2{i})*nonCO2_H2O;
                else
                    Inlet.Flow2.(block.Spec2{i}) = block.Flow2Spec.(block.Spec2{i})*nonCO2_H2O;
                end
            else Inlet.Flow2.(block.Spec2{i}) = 0;
            end
        end
end

% Fuel/Steam
for i = 1:1:length(block.Spec2)
    if isfield(block.Flow2Spec,block.Spec2{i})
        Inlet.Flow2.(block.Spec2{i}) = block.Flow2Spec.(block.Spec2{i})*block.Flow2_IC;%flow rate of every species entering the fuel/steam side (or reformer if there is one)
    else Inlet.Flow2.(block.Spec2{i}) = 0;
    end
end
end %Ends function InletFlow

function block = Set_IC(block,Flow1,Flow2,Flow3)
if NetFlow(Flow1.Inlet)>0
    Cp.flow1 = SpecHeat(Flow1.Inlet);
else
    Cp.flow1 = SpecHeat(Flow1.Outlet);
end
Cp.flow2 = SpecHeat(Flow2.Outlet);
if ~isempty(Flow3)
    Cp.flow3 = SpecHeat(Flow3.Outlet);
end
switch block.Reformer
    case 'internal'
        NumOfStates = (6 + length(block.Spec1) + 2*length(block.Spec2) + 1)*block.nodes + 2; % 6 temperatures, anode & cathode & reformer & current at each node and 2 states for anode/cathode pressure
  case {'direct';'external';'adiabatic';'none'}
        NumOfStates = (4 + length(block.Spec1) + length(block.Spec2) + 1)*block.nodes +2; % 4 temperatures, anode & cathode & current at each node and 2 states for anode/cathode pressure
end
block.IC = ones(NumOfStates,1); %
block.UpperBound = inf*ones(NumOfStates,1);
block.LowerBound = zeros(NumOfStates,1); %need to make this -inf for current states
block.tC = block.IC; % time constant for derivative dY
block.Scale = block.IC;
switch block.Reformer
    case 'internal'
        n = 6*block.nodes;
        block.tC(4*block.nodes+1:5*block.nodes) = (block.Mass_plate(2)*block.C_plate(2));
        block.tC(5*block.nodes+1:6*block.nodes) = (block.Vol_flow(3)*Cp.flow3*block.Flow2_Pinit./(block.Ru*block.T.Flow3));
    case {'direct';'external';'adiabatic';'none'}
        n = 4*block.nodes;
end
block.Scale = block.Tstates(1:n);%temperature (K)
block.tC(1:block.nodes) = (block.Mass_plate(1)*block.C_plate(1));
block.tC(1+block.nodes:2*block.nodes) = (block.Vol_flow(1)*Cp.flow1*block.Flow1_Pinit./(block.Ru*block.T.Flow1));
block.tC(2*block.nodes+1:3*block.nodes) = (block.Vol_Elec*block.Density_Elec*block.C_Elec);
block.tC(3*block.nodes+1:4*block.nodes) = (block.Vol_flow(2)*Cp.flow2*block.Flow2_Pinit./(block.Ru*block.T.Flow2));
block.tC(1:n) = block.tC(1:n)-diag(block.HTconv)-diag(block.HTcond); %this accounts for the change in HT as temperature of the control volume changes. The change in HT helps balance the energy equation more than the change in enthalpy leaving.

for i = 1:1:length(block.Spec1)
    block.tC(n+1:n+block.nodes) = (block.Vol_flow(1)*block.Flow1_Pinit)./(block.T.Flow1*block.Ru);  % air/ oxidant 
    if any(Flow1.Outlet.(block.Spec1{i})==0)
        block.IC(n+1:n+block.nodes) = Flow1.Outlet.(block.Spec1{i})./max(NetFlow(Flow1.Outlet));
        block.Scale(n+1:n+block.nodes) = max(NetFlow(Flow1.Outlet)); n = n+block.nodes; %cathode flows
    else
        block.Scale(n+1:n+block.nodes) = Flow1.Outlet.(block.Spec1{i}); n = n+block.nodes; %cathode flows
    end
end

for i = 1:1:length(block.Spec2)
    X = Flow2.Outlet.(block.Spec2{i})./NetFlow(Flow2.Outlet);%concentration
    block.tC(n+1:n+block.nodes) = (block.Vol_flow(2)*block.Flow2_Pinit)./(block.T.Flow2*block.Ru); %fuel/steam
    if all(X==0)
        block.IC(n+1:n+block.nodes) = 0;
        block.Scale(n+1:n+block.nodes) = NetFlow(Flow2.Outlet); %anode flow
    elseif any(X<.01) %concentration less than 1%
        block.IC(n+1:n+block.nodes) = Flow2.Outlet.(block.Spec2{i})/max(Flow2.Outlet.(block.Spec2{i}));
        block.Scale(n+1:n+block.nodes) = max(Flow2.Outlet.(block.Spec2{i})); %max flow of this species (inlet or outlet probably)
    else
        block.Scale(n+1:n+block.nodes) = Flow2.Outlet.(block.Spec2{i}); %individual species flow
    end   
    n = n+block.nodes;
end

switch block.Reformer
    case 'internal'
        for i = 1:1:length(block.Spec2)
            block.tC(n+1:n+block.nodes) = (block.Vol_flow(3)*block.Flow2_Pinit)./(block.T.Flow3*block.Ru); % reformer
            if any(Flow3.Outlet.(block.Spec2{i})==0)
                block.IC(n+1:n+block.nodes) = Flow3.Outlet.(block.Spec2{i})./NetFlow(Flow3.Outlet);
                block.Scale(n+1:n+block.nodes) = NetFlow(Flow3.Outlet); n = n+block.nodes; %anode flows
            else
                block.Scale(n+1:n+block.nodes) = Flow3.Outlet.(block.Spec2{i}); n = n+block.nodes; %reformer flows
            end
        end
end
block.tC(n+1:n+block.nodes) = block.nodes/100;%  %current changing for voltage balance
block.Scale(n+1:n+block.nodes) = (block.Current.H2+block.Current.CO); %current
block.LowerBound(n+1:n+block.nodes) = -inf; n = n+block.nodes; %current
block.tC(n+1) = (block.Vol_flow(1)*block.nodes*block.Cells); %pressure
block.tC(n+2) = (block.Vol_flow(2)*block.nodes*block.Cells);  %pressure
block.Scale(n+1) = block.Flow1_Pinit;%pressure
block.Scale(n+2) = block.Flow2_Pinit;%pressure
end %Ends function Set_IC

function block = loadGeometry(block,flow3)
%%--Geometric Varibles  %
block.t_plate = 0.003;                   % [m] thickness of seperator plate
block.t_plate_wall =0.005;                  % [m] Thickness of the channel wall of the Fuel Seperator Plate
block.H_platechannels = 0.002;                       % [m] height of flow 1 channel
block.W_platechannels = 0.005;                       % [m] width of channel                       
block.Nu_flow = 4;                           %Nusselt # for square channel aspect ratio =3 & uniform temp
block.t_plate_wall(2)=0.005;                     % [m] Thickness of the channel wall of the 2nd half of bi-polarPlate
block.H_platechannels(2) = 0.002;                      % [m] height of flow 2 channel
block.W_platechannels(2) = 0.005;                      % [m] width of flow 2 channel
block.Nu_flow(2) = 4;                           %Nusselt # for square channel aspect ratio =3 & uniform temp   
%%---Bi-polar plate-------%
block.Density_plate = 2000;                                % [kg/m3]     density 
block.C_plate = .600;                                        % [kJ/(kg K)] specific heat of fuel seperator plate
block.k_plate = 5;   %25                                    % [W/(m K)]   conductivity of Fuel Seperator Plate
%%----- Flow 1 (fuel or steam)-----------%
block.k_flow(1) = 259E-3;                                 % (W/m*K) Thermal Conductivity of 50%H2 & 50%H2O at 1000K
block.k_flow(2) = 67E-3;                                              % (W/m*K) Thermal Conductivity of air at 1000K
if flow3 %%----- / ReformerChannel---
    block.k_flow(3) = 112.52E-3;                                          % (W/m*K) Thermal Conductivity of 5%CO, 35%CO2 15%H2, 45%H20 at 900K
    %split bi-polar plate in half
    block.t_plate = block.t_plate/2;                   % [m] thickness of seperator plate
    block.t_plate(2) = block.t_plate(1)/2;                     % [m] Thickness of 2nd half of bi-polarPlate
    block.Density_plate(2) = 2000;                                % [kg/m3]     density of 2nd half of bi-polarPlate
    block.C_plate(2) = .600;                                        % [kJ/(kg K)] specific heat of 2nd half of bi-polarPlate
    block.k_plate(2) = 5;%25;                                % [W/(m K)]   conductivity of 2nd half of bi-polarPlate
    block.H_platechannels(3) = 0.005;                   % (m) height of flow 3 channel
    block.W_platechannels(3) = 0.003;                  % (m) width of flow 3 channel 
    block.t_plate_wall(3) = .001;               % (m)Thickness of flow 3 channel wall   
end
%%%---%%% end of user defined variables    

%Dimensions (calculated from user specifications)
block.A_Cell = block.L_Cell*block.W_Cell; %Cell Area
block.A_Node = block.A_Cell/block.nodes; %node Area
block.L_node = block.L_Cell/block.columns; %Node length in meters
block.W_node = block.W_Cell/block.rows;  %Node Width in meters
block.A_Node_Surf = block.A_Node;                                    % [m^2] Surface Area

for i = 1:1:length(block.H_platechannels)
    block.Dh_Flow(i) = 4*(block.H_platechannels(i)*block.W_platechannels(i))/(2*(block.H_platechannels(i)+block.W_platechannels(i))); %(m) hydraulic diameter of channel
    block.CH_Flow(i) = block.W_node/(block.W_platechannels(i)+block.t_plate_wall(i)); %Number of channels
end
%% --- Bi-polar plate-------%
block.A_plate_elecCond = block.t_plate_wall(1)*block.L_node*block.CH_Flow(1);   % [m2] Conduction area between the flow1 and the electrolyte
block.A_plate_heatCond = (block.H_platechannels(1)*block.t_plate_wall(1) + (block.W_platechannels(1)+block.t_plate_wall(1))*block.t_plate(1))*block.CH_Flow(1); %[m^2] conduction area between nodes
block.L_plate_heatCond = block.H_platechannels(1);                                     % [m] Length of conduction between flow1 and electrolyte
block.Mass_plate = (block.H_platechannels(1)*block.t_plate_wall(1) + (block.W_platechannels(1)+block.t_plate_wall(1))*block.t_plate(1))*block.CH_Flow(1)*block.L_node*block.Density_plate(1);
%% -----Flow 1(fuel or steam)-----------%
block.flow_crossArea(1) = block.H_platechannels(1)*block.W_platechannels(1)*block.CH_Flow(1);              % [m2] Crossectional Area of flow 1 entrance
block.h_flow(1) = block.Nu_flow(1)*block.k_flow(1)/block.Dh_Flow(1);                     % [W/m2/K]  Convection coefficient between the flow 1 and plate
block.A_flow_plate(1) = (2*block.H_platechannels(1) + block.W_platechannels(1))*block.L_node*block.CH_Flow(1);       % [m2]  Area in common between flow 1 and Plate for convection
block.A_flow_elec(1) = block.W_platechannels(1)*block.L_node*block.CH_Flow(1);                  % [m2]  Area in common between flow 1 and Electrolyte for convection
block.Vol_flow(1) = block.H_platechannels(1)*block.W_platechannels(1)*block.L_node*block.CH_Flow(1);               % [m3]  control volume flow 1
%% --------Electrolyte-------------%
block.A_Elec_Cond =  block.W_node*block.t_Elec;                   % [m2] Conduction surface area of electrolyte
block.A_Elec_Heat_Cond = block.W_node*block.t_Elec;                    % [m2] Conduction surface area of electrolyte
block.Vol_Elec = block.t_Elec*block.L_node*block.W_node;              % [m3] volume of electrolyte   
 
if flow3 %% ----2nd half of bi-polar plate-------%  
    block.Vol_flow(3) = block.H_platechannels(3)*block.W_platechannels(3)*block.CH_Flow(3)*block.L_node; % (m^3) Volume of Reformer Channel in cell
    block.flow_crossArea(3) = block.H_platechannels(3)*block.W_platechannels(3)*block.CH_Flow(3);     % (m^2) Reformer Channel Area per node
    block.h_flow(3) = block.Nu_flow(1)*block.k_flow(3)/block.Dh_Flow(3);                     % [W/m2/K]  Convection coefficient between the anode gas and the Fuel Seperator plate   
    
    block.A_plate_elecCond(2) = block.t_plate_wall(2)*block.L_node*block.CH_Flow(2);                % [m2] Conduction area between the fuel seperator plate and the electrolyte
    block.A_plate_heatCond(2) = (block.H_platechannels(2)*block.t_plate_wall(2) + (block.W_platechannels(2)+block.t_plate_wall(2))*block.t_plate(2))*block.CH_Flow(2); %[m^2] conduction area between nodes
    block.L_plate_heatCond(2) = block.H_platechannels(2);                                    % [m] Length of conduction between the fuel seperator plate and electrolyte
    block.Mass_plate(2) = (block.H_platechannels(2)*block.t_plate_wall(2) + (block.W_platechannels(2)+block.t_plate_wall(2))*block.t_plate(2))*block.L_node*block.CH_Flow(2)*block.Density_plate(2);
end
%% -------Flow 2 (Air or O2)---------%
block.flow_crossArea(2) = block.H_platechannels(2)*block.W_platechannels(2)*block.CH_Flow(2);       % [m2] Crossectional Area of flow 2
block.h_flow(2) = block.Nu_flow(2)*block.k_flow(2)/block.Dh_Flow(2);                 % [W/m2/K]  Convection coefficient between the flow 2 and plate
block.A_flow_plate(2) = (2*block.H_platechannels(2) + block.W_platechannels(2))*block.L_node*block.CH_Flow(2);    % [m2]  Area in common between flow 2 and plate for convection
block.A_flow_elec(2) = block.W_platechannels(2)*block.CH_Flow(2)*block.L_node;                 % [m2]  Area in common between flow 2 and Electrolyte for convection
block.Vol_flow(2) = block.H_platechannels(2)*block.W_platechannels(2)*block.CH_Flow(2)*block.L_node;            % [m3]  control volume flow2
end%ends function loadGeometry

function dY = DynamicTemps(t,Y,block,Flow1,Flow2,Flow3,Inlet)
dY = 0*Y;
nodes = block.nodes;

[h,hs] = enthalpy(Y(1+2*nodes:3*nodes),{'H2','H2O','O2','CO','CO2'});
Power = abs(block.Voltage)*(block.Current.H2+block.Current.CO)/1000; %node power in kW
Qreaction = block.Current.H2/(2000*block.F).*(h.H2+.5*h.O2-h.H2O) + block.Current.CO/(2000*block.F).*(h.CO+.5*h.O2-h.CO2);
Qgen = Qreaction-Power;%kW of heat generated by electrochemistry (per node & per cell)
switch block.FCtype%ion transport across membrane (total enthalpy)
    case 'SOFC'%ion crosses from flow 1 to flow 2 in fuel cell mode and from flow 2 to flow 1 in electrolyzer mode. Current is negative in electrolyzer mode
        Qion = (block.Current.H2+block.Current.CO)/(4000*block.F).*hs.O2; %O2 ion crossing over (kW)
    case 'MCFC'
        Qion = (block.Current.H2+block.Current.CO)*(1/(4000*block.F).*hs.O2 + 1/(2000*block.F).*hs.CO2);% O2 & CO2 ion crossing over
end

QT = block.HTconv*Y + block.HTcond*Y + block.HTrad*(Y.^4);

Flow1.Outlet.T = Y(nodes+1:2*nodes);
for j = 1:1:length(block.Flow1Dir(1,:));%1:columns
    k = block.Flow1Dir(:,j);
    if j~=1
        Flow1.Inlet.T(k,1) = Flow1.Outlet.T(kprev);
    end
    kprev = k;
end

Flow2.Outlet.T = Y(3*nodes+1:4*nodes);
for j = 1:1:length(block.Flow2Dir(1,:));%1:columns
    k = block.Flow2Dir(:,j);
    if j~=1
        Flow2.Inlet.T(k,1) = Flow2.Outlet.T(kprev);
    end
    kprev = k;
end

%energy flows & sepcific heats
Hout1 = enthalpy(Flow1.Outlet);
Hin1 = enthalpy(Flow1.Inlet);
Hout2 = enthalpy(Flow2.Outlet);
Hin2 = enthalpy(Flow2.Inlet);

if any(strcmp(block.Reformer,{'internal';'direct'})) && block.Recirc.Flow2>0 % Only during the first run with unhumidified fuel, find fuelmix temperature
    error2 = 1;
    Cp = SpecHeat(Flow2.Outlet); 
    Cp = Cp(end);
    Hmix = enthalpy(Inlet.Flow2) + sum(Hout2(block.Flow2Dir(:,end)))*block.Recirc.Flow2*block.Cells;
    netflow = NetFlow(Inlet.Mixed);
    while abs(error2)>1e-4
        error2 = (Hmix - enthalpy(Inlet.Mixed))./(Cp*netflow);                             %Adjusting the error in temperature based on known enthalpy and specific heat of the cold side
        Inlet.Mixed.T = Inlet.Mixed.T + .75*error2;                                   %Subtraction of a portion of the T_error from cold outlet temp to get closer to the actual temp
    end
end

switch block.Reformer
    case 'internal'
        Flow3.Outlet.T = Y(5*nodes+1:6*nodes);
        k = block.Flow3Dir(:,1);
        if block.Recirc.Flow2>0
            Flow3.Inlet.T(k,1) = Inlet.Mixed.T;
        else Flow3.Inlet.T(k,1) = Inlet.Flow2.T;
        end
        for j = 1:1:length(block.Flow3Dir(1,:));%1:columns
            k = block.Flow3Dir(:,j);
            if j~=1
                Flow3.Inlet.T(k,1) = Flow3.Outlet.T(kprev);
            end
            kprev = k;
        end
        Hin3 = enthalpy(Flow3.Inlet);
        Hout3 = enthalpy(Flow3.Outlet);
        
        k = block.Flow2Dir(:,1);
        k2 = block.Flow3Dir(:,end);
        Flow2.Inlet.T(k,1) = Flow3.Outlet.T(k2,1);
        Hin2 = enthalpy(Flow2.Inlet);
        if block.ClosedCathode %%energy balance
            Qimbalance = sum((Hin1(block.Flow1Dir(:,1))) - sum(Hout1(block.Flow1Dir(:,end)))) + sum(Hin3(block.Flow3Dir(:,1)))  - sum(Hout2(block.Flow2Dir(:,end))) - sum(Power);
            Power = Power + Qimbalance*Power./sum(Power);
            Qreaction = Power + Qgen;
        end
    case {'adiabatic';'direct';'external';'none';}
        k = block.Flow2Dir(:,1);
        Flow2.Inlet.T(k,1) = Inlet.Flow2.T;
        Hin2 = enthalpy(Flow2.Inlet);
        if block.ClosedCathode %%energy balance
            Qimbalance = sum((Hin1(block.Flow1Dir(:,1))) - sum(Hout1(block.Flow1Dir(:,end)))) + sum(Hin2(block.Flow2Dir(:,1)))  - sum(Hout2(block.Flow2Dir(:,end))) - sum(Power);
            Power = Power + Qimbalance*Power./sum(Power);
            Qreaction = Power + Qgen;
        end
end

dY(1:nodes)= QT(1:nodes)./block.tC(1:nodes);  %Bi-polar Plate
dY(1+nodes:2*nodes)= (QT(1+nodes:2*nodes) + Hin1 - Hout1 - Qion)./block.tC(1+nodes:2*nodes); %air/oxidant
dY(1+2*nodes:3*nodes)= (QT(1+2*nodes:3*nodes)+ Qgen)./block.tC(2*nodes+1:3*nodes); %Electrolyte Plate
dY(1+3*nodes:4*nodes)= (QT(1+3*nodes:4*nodes) + Hin2 - Hout2 + Qion - Qreaction)./block.tC(1+3*nodes:4*nodes);  %Fuel/steam
switch block.Reformer
    case 'internal'
        dY(1+4*nodes:5*nodes)= QT(1+4*nodes:5*nodes)./block.tC(4*nodes+1:5*nodes);  %split bi-polar plate
        dY(1+5*nodes:6*nodes)= (block.RefSpacing*QT(1+5*nodes:6*nodes) + Hin3 - Hout3)./block.tC(1+5*nodes:6*nodes);  %Fuel Reformer Channels
end
end %Ends function DynamicTemps