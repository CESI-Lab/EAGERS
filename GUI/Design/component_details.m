function component_details(varargin)
% This function updates the parameter options for the selected generator
% The user can then make changes to this particular generator and save
% those changes.
if nargin && ischar(varargin{1})
    f= str2func(varargin{1});
    f(varargin{2:length(varargin)})
else
    create_gen_details(varargin{1})
end

function create_gen_details(handles)
global testSystems SYSINDEX GENINDEX
set(handles.uipanelLibrary,'Visible','off');
set(handles.uipanelGenSpec,'Visible','on');
set(handles.Library,'Visible','on');
set(handles.saveSystem,'Visible','on');
set(handles.pushbuttonRemove,'Visible','on');
% if GENINDEX > 0
    gen = testSystems(SYSINDEX).Generator(GENINDEX);
% else
%     switch GENINDEX
%         case 0
%             Gen = struct('Type', 'None', ...
%                 'Name', 'None');
%         case -1
%             Gen = struct('Type', 'Heating Demands', ...
%                 'Name', 'Heating Demands', ...
%                 'Demand', 100);
%         case -2
%             Gen = struct('Type', 'Hot Water Demands', ...
%                 'Name', 'Hot Water Demands', ...
%                 'Demand', 100);
%         case -3
%             Gen = struct('Type', 'Cooling Demands', ...
%                 'Name', 'Cooling Demands', ...
%                 'Demand', 100);
%     end
% end
to_del=get(handles.uipanelGenSpec,'Children');
for i= 1:length(to_del)
    delete(to_del(i));
end
create_item(handles.uipanelGenSpec,'edit',gen.Name,[28 34.5 75 2],15,'bold',[],'CompName',[]);
switch gen.Type
    case 'Utility'
        if strcmp(gen.Source,'Electricity')            
            create_item(handles.uipanelGenSpec,'text','(1) Off-Peak Rate ',[47 32 25 2],12,'bold',[],'textEdit1',[]);
            create_item(handles.uipanelGenSpec,'text','Electric Charge ($/kWh) ',[57 30 32 1.75],11.5,'normal',[],'textEdit2',[]);
            create_item(handles.uipanelGenSpec,'text','Electric Charge ($/kWh) ',[57 30 32 1.75],11.5,'normal',[],'textEdit2',[]);
            create_item(handles.uipanelGenSpec,'text','Demand Charge ($/kW) ',[57 28 32 1.75],11.5,'normal',[],'textEdit3',[]);
            create_item(handles.uipanelGenSpec,'text','(2) Partial-Peak Rate ',[47 26 31 1.75],12,'bold',[],'textEdit4',[]);
            create_item(handles.uipanelGenSpec,'text','Electric Charge ($/kWh) ',[57 24 32 1.75],11.5,'normal',[],'textEdit5',[]);
            create_item(handles.uipanelGenSpec,'text','Demand Charge ($/kW) ',[57 22 32 1.75],11.5,'normal',[],'textEdit6',[]);
            create_item(handles.uipanelGenSpec,'text','(3) Peak Rate ',[47 20 20 1.75],12,'bold',[],'textEdit7',[]);
            create_item(handles.uipanelGenSpec,'text','Electric Charge ($/kWh) ',[57 18 32 1.75],11.5,'normal',[],'textEdit8',[]);
            create_item(handles.uipanelGenSpec,'text','Demand Charge ($/kW) ',[57 16 32 1.75],11.5,'normal',[],'textEdit9',[]);
            
            create_item(handles.uipanelGenSpec,'edit',num2str(gen.VariableStruct.SumRates(1,1)),[90 30 15 1.75],10,'normal',[],'compText1',[]);
            create_item(handles.uipanelGenSpec,'edit',num2str(gen.VariableStruct.SumRates(1,2)),[90 28 15 1.75],10,'normal',[],'compText2',[]);
            create_item(handles.uipanelGenSpec,'edit',num2str(gen.VariableStruct.SumRates(2,1)),[90 24 15 1.75],10,'normal',[],'compText3',[]);
            create_item(handles.uipanelGenSpec,'edit',num2str(gen.VariableStruct.SumRates(2,2)),[90 22 15 1.75],10,'normal',[],'compText4',[]);
            create_item(handles.uipanelGenSpec,'edit',num2str(gen.VariableStruct.SumRates(3,1)),[90 18 15 1.75],10,'normal',[],'compText5',[]);
            create_item(handles.uipanelGenSpec,'edit',num2str(gen.VariableStruct.SumRates(3,2)),[90 16 15 1.75],10,'normal',[],'compText6',[]);
            
            
            hrs = cell(24,1);
            for i = 1:1:length(hrs)
                hrs{i} = num2str(i);
            end
            colwidth = 21*ones(1,24);
            colwidth = num2cell(colwidth);
            
            create_table(handles,gen.VariableStruct.SumRateTable,hrs,colwidth,{'Sun';'Mon';'Tue';'Wed';'Thu';'Fri';'Sat';},[2 1 102 10.5],'uitableEffCurve',true(1,24));
            
            
            create_item(handles.uipanelGenSpec,'checkbox','Summer Rates',[1 33 30 1],12,'normal',1,'sum_rate','checkboxSeasonal');   
            create_item(handles.uipanelGenSpec,'checkbox','Winter Rates',[1 31.5 30 1],12,'normal',0,'win_rate','checkboxSeasonal');   
            create_item(handles.uipanelGenSpec,'text','Start Month',[1 29 20 1.75],11.5,'normal',[],'textEdit11',[]);
            create_item(handles.uipanelGenSpec,'text','Start Day',[22 29 20 1.75],11.5,'normal',[],'textEdit12',[]);
            create_item(handles.uipanelGenSpec,'edit',num2str(gen.VariableStruct.SumStartMonth),[4 26.5 14 1.75],10,'normal',[],'compText7',[]);
            create_item(handles.uipanelGenSpec,'edit',num2str(gen.VariableStruct.SumStartDay),[25 26.5 14 1.75],10,'normal',[],'compText8',[]);
            
            create_item(handles.uipanelGenSpec,'text','Grid Sell Back',[2 24.5 30 1.5],12,'bold',[],'textEdit10',[]);            
            
            if testSystems(SYSINDEX).Generator(GENINDEX).VariableStruct.SellBackRate == -1
                v1=0;v2=1;v3=0;
            elseif testSystems(SYSINDEX).Generator(GENINDEX).VariableStruct.MinImportThresh>=0 && testSystems(SYSINDEX).Generator(GENINDEX).VariableStruct.SellBackRate==0
                v1=1;v2=0;v3=0;
            else
                v1=0;v2=0;v3=1;
            end
            
            create_item(handles.uipanelGenSpec,'radiobutton','None',[2 23 10 1.5],10,'normal',v1,'None','radiobuttonGridsellback');
            create_item(handles.uipanelGenSpec,'radiobutton','% Purchase',[2 21 18 1.5],10,'normal',v2,'SellbackAsPerc','radiobuttonGridsellback');
            create_item(handles.uipanelGenSpec,'radiobutton','Fixed Rate',[2 19 18 1.5],10,'normal',v3,'FixedRate','radiobuttonGridsellback');

            if v1 == 1
                create_item(handles.uipanelGenSpec,'text','Minimum Import (kW)',[10 16 30 1.75],12,'normal',[],'textEdit11',[]);
            else
                create_item(handles.uipanelGenSpec,'text','Maximum Export (kW)',[10 16 30 1.75],12,'normal',[],'textEdit11',[]);
            end
            
            create_item(handles.uipanelGenSpec,'edit',num2str(abs(testSystems(SYSINDEX).Generator(GENINDEX).VariableStruct.MinImportThresh)),[20 14 12 1.5],10,'normal',[],'maxSellback',[]);
            create_item(handles.uipanelGenSpec,'edit',num2str(testSystems(SYSINDEX).Generator(GENINDEX).VariableStruct.SellBackPerc),[20 21 8 1.5],10,'normal',[],'editTariffs',[]);
            create_item(handles.uipanelGenSpec,'edit',num2str(max(0,testSystems(SYSINDEX).Generator(GENINDEX).VariableStruct.SellBackRate)),[20 19 8 1.5],10,'normal',[],'editSellbackRate',[]);
            
            create_item(handles.uipanelGenSpec,'text','$/kWh',[28 18.75 10 1.75],11.5,'normal',[],'textEdit12',[]);
        else % Gen.Source: NG (Natural Gas)
            create_item(handles.uipanelGenSpec,'text','Fuel Rate ($/MMBTU)',[60 32 30 1.75],12,'normal',[],'textEdit1',[]);
            create_item(handles.uipanelGenSpec,'edit',num2str(gen.VariableStruct.Rate(1)),[90 32 15 1.75],10,'normal',[],'compText1',[]);
        end
    case 'AC_DC'
        create_item(handles.uipanelGenSpec,'text','AC to DC (%)',[70 32 20 1.75],12,'normal',[],'textEdit1',[]);
        create_item(handles.uipanelGenSpec,'text','DC to AC (%)',[70 30 30 1.75],12,'normal',[],'textEdit2',[]);
        create_item(handles.uipanelGenSpec,'text','Capacity',[70 28 30 1.75],12,'normal',[],'textEdit3',[]);
        
        create_item(handles.uipanelGenSpec,'edit',num2str(gen.VariableStruct.AC_to_DC_eff),[90 32 15 1.75],10,'normal',[],'compText1',[]);
        create_item(handles.uipanelGenSpec,'edit',num2str(gen.VariableStruct.DC_to_AC_eff),[90 30 15 1.75],10,'normal',[],'compText2',[]);
        create_item(handles.uipanelGenSpec,'edit',num2str(gen.VariableStruct.Capacity),[90 28 15 1.75],10,'normal',[],'compText3',[]);
        
    case {'CHP Generator';'Electric Generator';'Heater';'Chiller';}
        if strcmp(gen.Type,'Heater')
            OutNames = {'Capacity';'Heat';};
            Data = [gen.Output.Capacity,gen.Output.Heat;];
            pos = [16 3 19.3 15.7143];
            LB = gen.VariableStruct.Startup.Heat(end);
        elseif strcmp(gen.Type,'Chiller')
            OutNames = {'Capacity';'Cooling';};
            pos = [10 3 20 16];
            Data = [gen.Output.Capacity,gen.Output.Cooling;];
            LB = gen.VariableStruct.Startup.Cooling(end);
        elseif strcmp(gen.Type,'Electric Generator')
            OutNames = {'Capacity';'Electricity';};
            pos = [10 3 20 16];
            if isfield(gen.Output,'Electricity')
                Data = [gen.Output.Capacity,gen.Output.Electricity;];
                LB = gen.VariableStruct.Startup.Electricity(end);
            elseif isfield(gen.Output,'DirectCurrent')
                Data = [gen.Output.Capacity,gen.Output.DirectCurrent;];
                LB = gen.VariableStruct.Startup.DirectCurrent(end);
            end
        elseif strcmp(gen.Type,'CHP Generator') 
            pos = [8 3 30 16];
            if isfield(gen.Output,'Electricity')
                Data = [gen.Output.Capacity,gen.Output.Electricity,gen.Output.Heat;];
                LB = gen.VariableStruct.Startup.Electricity(end);
            elseif isfield(gen.Output,'DirectCurrent')
                Data = [gen.Output.Capacity,gen.Output.DirectCurrent,gen.Output.Heat;];
                LB = gen.VariableStruct.Startup.DirectCurrent(end);
            end
            OutNames = {'Capacity';'Electricity';'Heat';};
        end
        create_table(handles,Data,OutNames,'auto',{},pos,'uitableEffCurve',true(1,length(Data(1,:))))
        
        create_item(handles.uipanelGenSpec,'text','Capacity (kW)',[70 32 20 1.75],12,'normal',[],'textEdit1',[]);
        create_item(handles.uipanelGenSpec,'text','Minimum Output (kW)',[60 30 30 1.75],12,'normal',[],'textEdit2',[]);
        create_item(handles.uipanelGenSpec,'text','Ramp Rate (kW / hr)',[60 28 30 1.75],12,'normal',[],'textEdit3',[]);
        create_item(handles.uipanelGenSpec,'text','Startup Cost ($/start)',[60 26 30 1.75],12,'normal',[],'textEdit4',[]);
        create_item(handles.uipanelGenSpec,'text','Nat. Freq.',[1 32.5 15 1.65],12,'normal',[],'textEdit5',[]);  
        create_item(handles.uipanelGenSpec,'text','Damping',[32 32.5 13 1.65],12,'normal',[],'textEdit6',[]);
        create_item(handles.uipanelGenSpec,'text','Energy Source',[70 23.5 22 1.75],12,'bold',[],'textEdit7',[]);
        
        create_item(handles.uipanelGenSpec,'edit',num2str(gen.Size),[90 32 15 1.75],10,'normal',[],'compText1',[]);
        create_item(handles.uipanelGenSpec,'edit',num2str(LB),[90 30 15 1.75],10,'normal',[],'compText2',[]);
        if isfield(gen.VariableStruct,'dX_dt')
            create_item(handles.uipanelGenSpec,'edit',num2str(gen.VariableStruct.dX_dt),[90 28 15 1.75],10,'normal',[],'compText3',[]);
            p = eig(gen.VariableStruct.StateSpace.A);
            w_0 = sqrt(real(p(1))^2 + imag(p(1))^2);
            zeta = -real(p(1)+p(2))/(2*w_0);
            create_item(handles.uipanelGenSpec,'edit',num2str(w_0),[16 32.5 14 1.75],10,'normal',[],'compText5',[]);
            create_item(handles.uipanelGenSpec,'edit',num2str(zeta),[45 32.5 14 1.75],10,'normal',[],'compText6',[]);
        else
            create_item(handles.uipanelGenSpec,'edit','--',[90 28 15 1.75],10,'normal',[],'compText3',[]);
            create_item(handles.uipanelGenSpec,'edit','--',[16 32.5 14 1.75],10,'normal',[],'compText5',[]);
            create_item(handles.uipanelGenSpec,'edit','--',[45 32.5 14 1.75],10,'normal',[],'compText6',[]);
        end
        if isfield(gen.VariableStruct,'StartCost')
            create_item(handles.uipanelGenSpec,'edit',num2str(gen.VariableStruct.StartCost),[90 26 15 1.75],10,'normal',[],'compText4',[]);
        else
            create_item(handles.uipanelGenSpec,'edit','--',[90 26 15 1.75],10,'normal',[],'compText4',[]);
        end
        
        if strcmp(gen.Source, 'NG')
            v1=1;v2=0;v3=0;v4=0;
        elseif strcmp(gen.Source, 'Oil')
            v1=0;v2=1;v3=0;v4=0;
        elseif strcmp(gen.Source, 'Electricity')
            v1=0;v2=0;v3=1;v4=0;
        elseif strcmp(gen.Source, 'Heat')
            v1=0;v2=0;v3=0;v4=1;
        end
        
        if strcmp(gen.Type,'Heater')
            create_item(handles.uipanelGenSpec,'radiobutton','Natural Gas',[60 22 25 1.75],10.5,'normal',v1,'NatGas','radioSource');
            create_item(handles.uipanelGenSpec,'radiobutton','Electricity',[85 22 18 1.75],10.5,'normal',v3,'Electricity','radioSource');
        elseif strcmp(gen.Type,'Chiller')
            create_item(handles.uipanelGenSpec,'radiobutton','Electricity',[60 22 25 1.75],10.5,'normal',v3,'Electricity','radioSource');
            create_item(handles.uipanelGenSpec,'radiobutton','Heat',[85 22 18 1.75],10.5,'normal',v4,'Heat','radioSource');
        else
            create_item(handles.uipanelGenSpec,'radiobutton','Natural Gas',[60 22 25 1.75],10.5,'normal',v1,'NatGas','radioSource')
            create_item(handles.uipanelGenSpec,'radiobutton','Diesel',[85 22 18 1.75],10.5,'normal',v2,'Diesel','radioSource')
        end
            
        %%Need to figure out plotting
        
        handles.ResponseRate = axes('Units','characters',...
            'Position', [9,22,40,10],'NextPlot','add',...
            'Tag', 'ResponseRate',...
            'Parent', handles.uipanelGenSpec,...
            'Visible','on');

        handles.EffCurve = axes('Units','characters',...
            'Position', [58,3.5,44,15],'NextPlot','add',...
            'Tag', 'EffCurve',...
            'Parent', handles.uipanelGenSpec,...
            'Visible','on');

        plotGenEfficiency(gen,handles)
        ss_response(gen,handles);

    case 'Solar'
        create_item(handles.uipanelGenSpec,'text','Location',[66 32 14 1.75],12,'normal',[],'textEdit1',[]);
        create_item(handles.uipanelGenSpec,'text','Size (kW)',[75 28 14 1.75],12,'normal',[],'textEdit2',[]);
        create_item(handles.uipanelGenSpec,'text','Size (m^2)',[73.5 26 16 1.75],12,'normal',[],'textEdit3',[]);
        create_item(handles.uipanelGenSpec,'text','Conversion Efficiency (%)',[55 22 34 1.75],11.5,'normal',[],'textEdit4',[]);
        create_item(handles.uipanelGenSpec,'text','Azimuth Angle (Degrees) ',[57 20 32 1.75],11.5,'normal',[],'textEdit5',[]);
        create_item(handles.uipanelGenSpec,'text','Tilt Angle (Degrees)',[63 18 26 1.75],11.5,'normal',[],'textEdit6',[]);
        create_item(handles.uipanelGenSpec,'text','Solar Type',[15 28.1 17 1.75],12,'bold',[],'textEdit7',[]);
        create_item(handles.uipanelGenSpec,'text','Solar Tracking',[12 24.1 22 1.75],12,'bold',[],'textEdit8',[]);

        create_item(handles.uipanelGenSpec,'edit',num2str(gen.VariableStruct.Size),[90 28 15 1.75],10,'normal',[],'compText1',[]);
        create_item(handles.uipanelGenSpec,'edit',num2str(gen.VariableStruct.Sizem2),[90 26 15 1.75],10,'normal',[],'compText2',[]);
        create_item(handles.uipanelGenSpec,'edit',num2str(gen.VariableStruct.Eff),[90 22 15 1.75],10,'normal',[],'compText3',[]);
        create_item(handles.uipanelGenSpec,'edit',num2str(gen.VariableStruct.Azimuth),[90 20 15 1.75],10,'normal',[],'compText4',[]);
        create_item(handles.uipanelGenSpec,'edit',num2str(gen.VariableStruct.Tilt),[90 18 15 1.75],10,'normal',[],'compText5',[]);

        stateName = {'Alabama';'Alaska';'Arizona';'Arkansas';'California';'Colorado';'Connecticut';'Delaware';'Florida';'Georgia';'Hawaii';'Idaho';'Illinois';'Indiana';'Iowa';'Kansas';
         'Kentucky';'Louisiana';'Maine';'Maryland';'Massachusetts';'Michigan';'Minnesota';'Mississippi';'Missouri';'Montana';'Nebraska';'Nevada';'New Hampshire';'New Jersey';
         'New Mexico';'New York';'North Carolina';'North Dakota';'Ohio';'Oklahoma';'Oregon';'Pennsylvania';'Rhode Island';'South Carolina';'South Dakota';'Tennessee';'Texas';
         'Utah';'Vermont';'Virginia';'Washington';'West Virginia';'Wisconsin';'Wyoming';};
        stateNum = find(strcmp(gen.VariableStruct.State,stateName));
        create_item(handles.uipanelGenSpec,'popup',stateName,[80 32 25 1.75],10,'normal',stateNum,'popupSolar',[]);
        if strcmp(gen.VariableStruct.PVtype,'flat')
            val1 = 1;
            val2 = 0;
        else
            val1 = 0;
            val2 = 1;
        end
        create_item(handles.uipanelGenSpec,'radiobutton','Flat Panel',[5 27 16 1.5],11,'normal',val1,'Flat','radiobuttonSType');
        create_item(handles.uipanelGenSpec,'radiobutton','Concentrated',[23 27 20 1.5],11,'normal',val2,'Concentrated','radiobuttonSType');

        if strcmp(gen.VariableStruct.Tracking,'fixed')
            val1 = 1;
            val2 = 0;
            val3 = 0;
        elseif strcmp(gen.VariableStruct.Tracking,'1axis')
            val1 = 0;
            val2 = 1;
            val3 = 0;
        else
            val1 = 0;
            val2 = 0;
            val3 = 1;
        end
        create_item(handles.uipanelGenSpec,'radiobutton','Fixed',[1 23 10 1.5],11,'normal',val1,'Fixed','radiobuttonSTracking');
        create_item(handles.uipanelGenSpec,'radiobutton','Single Axis',[13 23 17 1.5],11,'normal',val2,'SingleAxis','radiobuttonSTracking');            
        create_item(handles.uipanelGenSpec,'radiobutton','Dual Axis',[32 23 15 1.5],11,'normal',val3,'DualAxis','radiobuttonSTracking');        

        row = {'DC rating';'Inverter/Transformer';'Mismatch';'Diodes/Connection';'DC wiring';...
               'AC wiring';'Soiling';'System availability';'Shading';'Sun-Tracking';'Age';};
        create_table(handles,gen.VariableStruct.Data,{'Value';'Min';'Max'},'auto',row,[20 1.5 71.13 15.72],'uitableDCAC',logical([1 0 0]));

    case 'Electric Storage'
        create_item(handles.uipanelGenSpec,'text','Size (kWh)',[70 32 20 1.75],12,'normal',[],'textEdit1',[]);
        create_item(handles.uipanelGenSpec,'text','Voltage (V)',[70 30 20 1.75],12,'normal',[],'textEdit2',[]);
        create_item(handles.uipanelGenSpec,'text','Max Depth of Discharge (%)',[69 27 20 3],12,'normal',[],'textEdit3',[]);
        create_item(handles.uipanelGenSpec,'text','Charging Internal Resistance (mOhms @ 100A)',[64 22 25 4],12,'normal',[],'textEdit4',[]);
        create_item(handles.uipanelGenSpec,'text','Discharging Internal Resistance  (mOhms @ 100A)',[61.5 16.5 28 4],12,'normal',[],'textEdit5',[]);
        create_item(handles.uipanelGenSpec,'text','Peak Charge Rate (C)',[9 32 31 1.75],12,'normal',[],'textEdit6',[]);
        create_item(handles.uipanelGenSpec,'text','Peak Disharge Rate (C)',[6 30 35 1.75],12,'normal',[],'textEdit7',[]);
        create_item(handles.uipanelGenSpec,'text','Self Discharge Rate (%)',[6.5 28 34 1.75],12,'normal',[],'textEdit8',[]);

        create_item(handles.uipanelGenSpec,'edit',num2str(gen.Size),[90 32 15 1.75],10,'normal',[],'compText1',[]);
        create_item(handles.uipanelGenSpec,'edit',num2str(gen.VariableStruct.Voltage),[90 30 15 1.75],10,'normal',[],'compText2',[]);
        create_item(handles.uipanelGenSpec,'edit',num2str(gen.VariableStruct.MaxDOD),[90 27.5 15 1.75],10,'normal',[],'compText3',[]);
        create_item(handles.uipanelGenSpec,'edit',num2str(gen.VariableStruct.ChargeResist),[90 23 15 1.75],10,'normal',[],'compText4',[]);
        create_item(handles.uipanelGenSpec,'edit',num2str(gen.VariableStruct.DischResist),[90 17.5 15 1.75],10,'normal',[],'compText5',[]);
        create_item(handles.uipanelGenSpec,'edit',num2str(gen.VariableStruct.PeakCharge),[41 32 15 1.75],10,'normal',[],'compText6',[]);
        create_item(handles.uipanelGenSpec,'edit',num2str(gen.VariableStruct.PeakDisch),[41 30 15 1.75],10,'normal',[],'compText7',[]);
        create_item(handles.uipanelGenSpec,'edit',num2str((gen.VariableStruct.SelfDischarge)*(31*24*100)),[41 28 15 1.75],10,'normal',[],'compText8',[]);

        create_table(handles,gen.VariableStruct.VoltCurve,{'State of Charge';'Voltage'},'auto','numbered',[17 10 32.13 14.5],'uitableBat',true(1,2));

    case 'Thermal Storage'
        create_item(handles.uipanelGenSpec,'text','Size (kWh)',[70 32 20 1.75],12,'normal',[],'textEdit1',[]);
        create_item(handles.uipanelGenSpec,'text','Size (L)',[72 30 20 1.75],12,'normal',[],'textEdit2',[]);
        create_item(handles.uipanelGenSpec,'text','T Hot(C)',[72 28 20 1.75],12,'normal',[],'textEdit3',[]);
        create_item(handles.uipanelGenSpec,'text','T Cold(C)',[71.5 26 20 1.75],12,'normal',[],'textEdit4',[]);
        create_item(handles.uipanelGenSpec,'text','Ramp Rate',[71 24 20 1.75],12,'normal',[],'textEdit5',[]);
        create_item(handles.uipanelGenSpec,'text','Charging Efficiency (%)',[10.75 32 32 1.75],12,'normal',[],'textEdit6',[]);
        create_item(handles.uipanelGenSpec,'text','Discharging Efficiency (%)',[6 30 38 1.75],12,'normal',[],'textEdit7',[]);
        create_item(handles.uipanelGenSpec,'text','Self Discharge(% lost per day)',[1.5 28 42 1.75],12,'normal',[],'textEdit8',[]);
        create_item(handles.uipanelGenSpec,'text','Fill Rate (L/min)',[15 26 34 1.75],12,'normal',[],'textEdit9',[]);
        create_item(handles.uipanelGenSpec,'text','Discharge Rate (L/min)',[10.25 24 34 1.75],12,'normal',[],'textEdit10',[]);
        
        create_item(handles.uipanelGenSpec,'edit',num2str(gen.Size),[90 32 15 1.75],10,'normal',[],'compText1',[]);
        create_item(handles.uipanelGenSpec,'edit',num2str(gen.VariableStruct.SizeLiter),[90 30 15 1.75],10,'normal',[],'compText2',[]);
        create_item(handles.uipanelGenSpec,'edit',num2str(gen.VariableStruct.Tcold),[90 28 15 1.75],10,'normal',[],'compText3',[]);
        create_item(handles.uipanelGenSpec,'edit',num2str(gen.VariableStruct.Thot),[90 26 15 1.75],10,'normal',[],'compText4',[]);
        create_item(handles.uipanelGenSpec,'edit',num2str(gen.VariableStruct.dX_dt),[90 24 15 1.75],10,'normal',[],'compText5',[]);
        create_item(handles.uipanelGenSpec,'edit',num2str(gen.VariableStruct.ChargeEff),[45 32 15 1.75],10,'normal',[],'compText6',[]);
        create_item(handles.uipanelGenSpec,'edit',num2str(gen.VariableStruct.DischargeEff),[45 30 15 1.75],10,'normal',[],'compText7',[]);
        create_item(handles.uipanelGenSpec,'edit',num2str(gen.VariableStruct.SelfDischarge),[45 28 15 1.75],10,'normal',[],'compText8',[]);
        create_item(handles.uipanelGenSpec,'edit',num2str(gen.VariableStruct.FillRate),[45 26 15 1.75],10,'normal',[],'compText9',[]);
        create_item(handles.uipanelGenSpec,'edit',num2str(gen.VariableStruct.DischRate),[45 24 15 1.75],10,'normal',[],'compText10',[]);
        
%     case 'Boiler'
%         OutNames = {'Capacity';'Steam';};
%         Data(:,1) = Gen.Output.Capacity;
%         Data(:,2) = Gen.Output.Steam;
%         editable = [true true];
%         createTable(handles,Data,OutNames,'auto',{},[16 3 19.3 15.7143],'uitableEffCurve',editable)
%         
%         create_item(handles.uipanelGenSpec,'text','Capacity (kW)',[70 32 20 1.75],12,'normal',[],'textEdit1')
%         
%         create_item(handles.uipanelGenSpec,'edit',num2str(Gen.Size),[90 32 15 1.75],10,'normal',[],'compText1')
%         
%         plotGenEfficiency(Gen,handles)
        
    case 'Heating Demands'

    case 'Hot Water Demands'

    case 'Cooling Demands'

    case 'AC/DC Conversion'

    case 'None'
        
end

function create_item(h,style,string,position,size,weight,value,tag,call)
quote = '''';
if isempty(call)
    call = eval(strcat('@(hObject,eventdata)component_edit(',quote,tag,quote,',hObject,eventdata,guidata(hObject))'));
else
    call = eval(strcat('@(hObject,eventdata)component_edit(',quote,call,quote,',hObject,eventdata,guidata(hObject))'));
end
uicontrol('Style',style,'String',string,'Parent',h,'Units','Characters',...
      'Position',position,'Value',value,'FontSize',size,'FontWeight',weight,...
      'Tag',tag,'Callback',call);

function create_table(handles,Data,ColName,ColWidth,RowName,position,tag,edit)
uitable('Parent',handles.uipanelGenSpec,'Units','characters','Position',position,...
     'FontSize',8,'Data',Data,'RowName',RowName,'ColumnName',ColName,...
     'ColumnWidth',ColWidth,'Tag',tag,'ColumnEditable',edit,...
     'CellEditCallback',strcat('component_edit(',tag,',hObject,eventdata,guidata(hObject))'));