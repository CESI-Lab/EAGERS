function component_edit(varargin)
global testSystems SYSINDEX GENINDEX
gen = testSystems(SYSINDEX).Generator(GENINDEX);
if ischar(varargin{1})
    f= str2func(varargin{1});
    testSystems(SYSINDEX).Generator(GENINDEX) = f(gen,varargin{2:end});
else
    f = str2func(get(varargin{1},'Call'));
    testSystems(SYSINDEX).Generator(GENINDEX) = f(gen,varargin);
end
testSystems(SYSINDEX).Design = [];%empty design day solution
end%ends function component_edit

function gen = compText1(gen, hObject, eventdata, handles)
if strcmp(gen.Type,'Utility')
    gen.VariableStruct.SumRates(1,1) = str2double(get(hObject,'String'));
elseif strcmp(gen.Type,'AC_DC')
    gen.VariableStruct.AC_to_DC_eff = str2double(get(hObject,'String'));
elseif ismember(gen.Type,{'CHP Generator';'Electric Generator';'Heater';'Chiller'}) 
    handles = guihandles;
    gen = update_component_spec(gen,'UB',str2double(get(hObject,'String')));
    p = eig(gen.VariableStruct.StateSpace.A);
    w_0 = sqrt(real(p(1))^2 + imag(p(1))^2);
    zeta = -real(p(1)+p(2))/(2*w_0);
    set(handles.compText5,'String',num2str(w0));
    set(handles.compText6,'String',num2str(zeta));
    ss_response(gen,handles);
elseif ismember(gen.Type,{'Electric Storage';'Thermal Storage';}) 
    gen.Size = str2double(get(hObject,'String'));
elseif strcmp(gen.Type,'Solar')
    gen.VariableStruct.Size = str2double(get(hObject,'String'));
end
end%Ends function compText1

function gen = compText2(gen, hObject, eventdata, handles)
if strcmp(gen.Type,'Utility')
    gen.VariableStruct.SumRates(1,2) = str2double(get(hObject,'String'));
elseif strcmp(gen.Type,'AC_DC')
    gen.VariableStruct.DC_to_AC_eff = str2double(get(hObject,'String'));
elseif ismember(gen.Type,{'CHP Generator';'Electric Generator';'Heater';'Chiller'}) 
    handles = guihandles;
    gen = update_component_spec(gen,'LB',str2double(get(hObject,'String')));
    p = eig(gen.VariableStruct.StateSpace.A);
    w0 = sqrt(real(p(1))^2 + imag(p(1))^2);
    zeta = -real(p(1)+p(2))/(2*w0);
    set(handles.compText5,'String',num2str(w0));
    set(handles.compText6,'String',num2str(zeta));
    ss_response(gen,handles);
elseif strcmp(gen.Type,'Solar')
    gen.VariableStruct.Sizem2 = str2double(get(hObject,'String'));
elseif strcmp(gen.Type,'Electric Storage')
    gen.VariableStruct.Voltage = str2double(get(hObject,'String'));
elseif strcmp(gen.Type,'Thermal Storage')
    gen.VariableStruct.SizeLiter = str2double(get(hObject,'String'));
end
end%Ends function compText2

function gen = compText3(gen, hObject, eventdata, handles)
if strcmp(gen.Type,'Utility')
    gen.VariableStruct.SumRates(2,1) = str2double(get(hObject,'String'));
elseif strcmp(gen.Type,'AC_DC')
    gen.VariableStruct.Capacity = str2double(get(hObject,'String'));
elseif ismember(gen.Type,{'CHP Generator';'Electric Generator';'Heater';'Chiller'})
    handles = guihandles;
    gen = update_component_spec(gen,'dX_dt',str2double(get(hObject,'String')));
    p = eig(gen.VariableStruct.StateSpace.A);
    w0 = sqrt(real(p(1))^2 + imag(p(1))^2);
    zeta = -real(p(1)+p(2))/(2*w0);
    set(handles.compText5,'String',num2str(w0));
    set(handles.compText6,'String',num2str(zeta));
    ss_response(gen,handles);
elseif strcmp(gen.Type,'Solar')
    gen.VariableStruct.Eff = str2double(get(hObject,'String'));
elseif strcmp(gen.Type,'Electric Storage')
    gen.VariableStruct.MaxDOD = str2double(get(hObject,'String'));
elseif strcmp(gen.Type,'Thermal Storage')
    gen.VariableStruct.Tcold = str2double(get(hObject,'String'));
end
end%Ends function compText3

function gen = compText4(gen, hObject, eventdata, handles)
if strcmp(gen.Type,'Utility')
    gen.VariableStruct.SumRates(2,2) = str2double(get(hObject,'String'));
elseif ismember(gen.Type,{'CHP Generator';'Electric Generator';'Heater';'Chiller';})
    gen.VariableStruct.StartCost = str2double(get(hObject,'String'));
elseif strcmp(gen.Type,'Solar')
    gen.VariableStruct.Azimuth = str2double(get(hObject,'String'));
elseif strcmp(gen.Type,'Electric Storage')
    gen.VariableStruct.ChargeResist = str2double(get(hObject,'String'));
elseif strcmp(gen.Type,'Thermal Storage')
    gen.VariableStruct.Thot = str2double(get(hObject,'String'));
end
end%Ends function compText4

function gen = compText5(gen, hObject, eventdata, handles)
if strcmp(gen.Type,'Utility')
    gen.VariableStruct.SumRates(3,1) = str2double(get(hObject,'String'));
elseif ismember(gen.Type,{'CHP Generator';'Electric Generator';'Heater';'Chiller';})
    handles = guihandles;
    gen = update_component_spec(gen,'w0',str2double(get(hObject,'String')));
    set(handles.compText3,'String',num2str(gen.VariableStruct.dX_dt));
    ss_response(gen,handles);
elseif strcmp(gen.Type,'Solar')
    gen.VariableStruct.Tilt = str2double(get(hObject,'String'));
elseif strcmp(gen.Type,'Electric Storage')
    gen.VariableStruct.DischResist = str2double(get(hObject,'String'));
elseif strcmp(gen.Type,'Thermal Storage')
    gen.VariableStruct.dX_dt = str2double(get(hObject,'String'));
end
end%Ends function compText5

function gen = compText6(gen, hObject, eventdata, handles)
if strcmp(gen.Type,'Utility')
    gen.VariableStruct.SumRates(3,2) = str2double(get(hObject,'String'));
elseif ismember(gen.Type,{'CHP Generator';'Electric Generator';'Heater';'Chiller';})
    handles = guihandles;
    gen = update_component_spec(gen,'zeta',str2double(get(hObject,'String')));
    p = eig(gen.VariableStruct.StateSpace.A);
    w0 = sqrt(real(p(1))^2 + imag(p(1))^2);
    zeta = -real(p(1)+p(2))/(2*w0);
    set(handles.compText5,'String',num2str(w0));
    ss_response(gen,handles);
elseif strcmp(gen.Type,'Electric Storage')
    gen.VariableStruct.PeakCharge = str2double(get(hObject,'String'));
elseif strcmp(gen.Type,'Thermal Storage')
    gen.VariableStruct.ChargeEff = str2double(get(hObject,'String'));
end
end%Ends function compText6

function gen = compText7(gen, hObject, eventdata, handles)
if strcmp(gen.Type,'Utility')
    if get(handles.sum_rate,'Value') == 1
        gen.VariableStruct.SumStartMonth = str2double(get(hObject,'String'));
    else
        gen.VariableStruct.WinStartMonth = str2double(get(hObject,'String'));
    end
elseif strcmp(gen.Type,'Electric Storage')
    gen.VariableStruct.PeakDisch = str2double(get(hObject,'String'));
elseif strcmp(gen.Type,'Thermal Storage')
    gen.VariableStruct.DischargeEff = str2double(get(hObject,'String'));
end
end%Ends function compText7

function gen = compText8(gen, hObject, eventdata, handles)
if strcmp(gen.Type,'Utility')
    if get(handles.sum_rate,'Value') == 1
        gen.VariableStruct.SumStartDay = str2double(get(hObject,'String'));
    else
        gen.VariableStruct.WinStartDay = str2double(get(hObject,'String'));
    end
elseif strcmp(gen.Type,'Electric Storage')
    gen.VariableStruct.SelfDischarge = str2double(get(hObject,'String'))/(31*24*100);
elseif strcmp(gen.Type,'Thermal Storage')
    gen.VariableStruct.SelfDischarge = str2double(get(hObject,'String'));
end
end%Ends function compText8

function gen = compText9(gen, hObject, eventdata, handles)
if strcmp(gen.Type,'Thermal Storage')
    gen.VariableStruct.FillRate = str2double(get(hObject,'String'));
end
end%Ends function compText9

function gen = compText10(gen, hObject, eventdata, handles)
if strcmp(gen.Type,'Thermal Storage')
    gen.VariableStruct.DischRate = str2double(get(hObject,'String'));
end
end%Ends function compText10

function gen = uitableEffCurve(gen,hObject, eventdata, handles)
Outputs = fieldnames(gen.Output);
handlesC=get(handles.uipanelGenSpec,'Children');
nOutput = eventdata.Indices;
newValue = eventdata.NewData;
if length(Outputs) == 5
    for i= 1:length(handlesC)
        if strcmp(get(handlesC(i),'Tag'),'checkboxSeasonal')
            if get(handlesC(i),'Value') == 1
                type = get(handlesC(i-1),'Value');
                if type == 1
                    gen.VariableStruct.SumRateTable(nOutput(1),nOutput(2)) = newValue; 
                else
                    gen.VariableStruct.WinRateTable(nOutput(1),nOutput(2)) = newValue;
                end
            else
                gen.VariableStruct.SumRateTable(nOutput(1),nOutput(2)) = newValue;
            end
        end
    end
else
    gen.Output.(Outputs{nOutput(2)})(nOutput(1)) = newValue;
end
end%Ends function uitableEffCurve

function gen = popupRates(gen,hObject,eventdata,handles)
m_d=[1 31; 2 28; 3 31; 4 30; 5 31; 6 30; 7 31; 8 31; 9 30; 10 31; 11 30; 12 31];
handlesC = get(handles.uipanelGenSpec,'Children');
if hObject.Value == 1
    for i=1:length(handlesC)
        if strcmp(get(handlesC(i),'Tag'),'compText1')
            set(handlesC(i),'String',num2str(gen.VariableStruct.SumRates(1,1)));
        elseif strcmp(get(handlesC(i),'Tag'),'compText2')
            set(handlesC(i),'String',num2str(gen.VariableStruct.SumRates(1,2)));
        elseif strcmp(get(handlesC(i),'Tag'),'compText3')
            set(handlesC(i),'String',num2str(gen.VariableStruct.SumRates(2,1)));
        elseif strcmp(get(handlesC(i),'Tag'),'compText4')
            set(handlesC(i),'String',num2str(gen.VariableStruct.SumRates(2,2)));
        elseif strcmp(get(handlesC(i),'Tag'),'compText5')
            set(handlesC(i),'String',num2str(gen.VariableStruct.SumRates(3,1)));
        elseif strcmp(get(handlesC(i),'Tag'),'compText6')
            set(handlesC(i),'String',num2str(gen.VariableStruct.SumRates(3,2)));
        elseif strcmp(get(handlesC(i),'Tag'),'uitableEffCurve')
            set(handlesC(i),'Data',gen.VariableStruct.SumRateTable);
        end
    end
    Start1 = gen.VariableStruct.SumStartMonth;
    Start2 = gen.VariableStruct.SumStartDay;
    if gen.VariableStruct.WinStartDay == 1
        End1 = gen.VariableStruct.WinStartMonth-1;
        End2 = m_d(End1,2);
    else
        End1 = gen.VariableStruct.WinStartMonth;
        End2 = gen.VariableStruct.WinStartDay - 1;
    end
    tag1 = 'popupDates1S';
    tag2 = 'popupDates2S';
    tag3 = 'popupDates3S';
    tag4 = 'popupDates4S';
else
    for i=1:length(handlesC)
        if strcmp(get(handlesC(i),'Tag'),'compText1')
            set(handlesC(i),'String',num2str(gen.VariableStruct.WinRates(1,1)));
        elseif strcmp(get(handlesC(i),'Tag'),'compText2')
            set(handlesC(i),'String',num2str(gen.VariableStruct.WinRates(1,2)));
        elseif strcmp(get(handlesC(i),'Tag'),'compText3')
            set(handlesC(i),'String',num2str(gen.VariableStruct.WinRates(2,1)));
        elseif strcmp(get(handlesC(i),'Tag'),'compText4')
            set(handlesC(i),'String',num2str(gen.VariableStruct.WinRates(2,2)));
        elseif strcmp(get(handlesC(i),'Tag'),'compText5')
            set(handlesC(i),'String',num2str(gen.VariableStruct.WinRates(3,1)));
        elseif strcmp(get(handlesC(i),'Tag'),'compText6')
            set(handlesC(i),'String',num2str(gen.VariableStruct.WinRates(3,2)));
        elseif strcmp(get(handlesC(i),'Tag'),'uitableEffCurve')
            set(handlesC(i),'Data',gen.VariableStruct.WinRateTable);
        end
    end
    Start1 = gen.VariableStruct.WinStartMonth;
    Start2 = gen.VariableStruct.WinStartDay;
    if gen.VariableStruct.SumStartDay == 1
        End1 = gen.VariableStruct.SumStartMonth-1;
        End2 = m_d(End1,2);
    else
        End1 = gen.VariableStruct.SumStartMonth;
        End2 = gen.VariableStruct.SumStartDay - 1;
    end
    tag1 = 'popupDates1W';
    tag2 = 'popupDates2W';
    tag3 = 'popupDates3W';
    tag4 = 'popupDates4W';
end
days=(1:m_d(Start1,2))';
component_details('create_popup',handles,{'1','2','3','4','5','6','7','8','9','10','11','12'},...
            [10 29 7 1],12,'normal',tag1,Start1,'popupDates')
component_details('create_popup',handles,days,[18 29 7 1],12,'normal',tag2,Start2,'popupDates')
clear days
days=(1:m_d(End1,2))';
component_details('create_popup',handles,{'1','2','3','4','5','6','7','8','9','10','11','12'},...
            [30 29 7 1],12,'normal',tag3,End1,'popupDates')
component_details('create_popup',handles,days,[38 29 7 1],12,'normal',tag4,End2,'popupDates')
end%Ends function popupRates

function gen = popupDates(gen,hObject,eventdata,handles)
handlesc = get(handles.uipanelGenSpec,'Children');
xW=[];
xS=[];

for i=1:(length(handlesc))
    if ~isempty(strfind(handlesc(i).Tag,'popupDates')) && ~isempty(strfind(handlesc(i).Tag,'S'))
        xS(:,end+1)=[get(handlesc(i),'Value');i];
    elseif ~isempty(strfind(handlesc(i).Tag,'popupDates')) && ~isempty(strfind(handlesc(i).Tag,'W'))
        xW(:,end+1)=[get(handlesc(i),'Value');i];
    end
end
Date = get(hObject,'Value');
Type = get(hObject,'Tag');
m_d=[1 31; 2 28; 3 31; 4 30; 5 31; 6 30; 7 31; 8 31; 9 30; 10 31; 11 30; 12 31];


if ismember(Type,{'popupDates1S','popupDates2S','popupDates3S','popupDates4S'})
    xS = xS(:,1:4);
    xS = fliplr(xS);
    if strcmp(Type,'popupDates1S')
        gen.VariableStruct.SumStartMonth = Date;
    elseif strcmp(Type,'popupDates2S')
        gen.VariableStruct.SumStartDay = Date;
    elseif strcmp(Type,'popupDates3S')
        if xS(1,4) >= m_d(xS(1,3),2)
            gen.VariableStruct.WinStartMonth = Date + 1;
            gen.VariableStruct.WinStartDay = 1;
            if gen.VariableStruct.WinStartMonth == 13
                gen.VariableStruct.WinStartMonth = 1;
            end
        else
            gen.VariableStruct.WinStartMonth = Date;
            gen.VariableStruct.WinStartDay = xS(1,4)+1;
        end
    elseif strcmp(Type,'popupDates4S')
        if m_d(xS(1,3),2) == Date
            gen.VariableStruct.WinStartMonth = xS(1,3) + 1;
            gen.VariableStruct.WinStartDay = 1;
            if gen.VariableStruct.WinStartMonth == 13
                gen.VariableStruct.WinStartMonth = 1;
            end
        else
            gen.VariableStruct.WinStartMonth = xS(1,3);
            gen.VariableStruct.WinStartDay = Date + 1;
        end
    end
    days=(1:m_d(xS(1,1),2))';
    if get(handlesc(xS(2,2)),'Value')>length(days)
        set(handlesc(xS(2,2)),'Value',days(end))
    end
    set(handlesc(xS(2,2)),'String',days);
    clear days
    days=(1:m_d(xS(1,3),2))';
    if get(handlesc(xS(2,4)),'Value')>length(days)
        set(handlesc(xS(2,4)),'Value',days(end))
    end
    set(handlesc(xS(2,4)),'String',days);
elseif ismember(Type,{'popupDates1W','popupDates2W','popupDates3W','popupDates4W'})
    xW = xW(:,1:4);
    xW = fliplr(xW);
    if strcmp(Type,'popupDates1W')
        gen.VariableStruct.WinStartMonth = Date;
    elseif strcmp(Type,'popupDates2W')
        gen.VariableStruct.WinStartDay = Date;
    elseif strcmp(Type,'popupDates3W')
        if xW(1,4) >= m_d(xW(1,3),2)
            gen.VariableStruct.SumStartMonth = Date + 1;
            gen.VariableStruct.SumStartDay = 1;
            if gen.VariableStruct.SumStartMonth == 13
                gen.VariableStruct.SumStartMonth = 1;
            end
        else
            gen.VariableStruct.SumStartMonth = Date;
            gen.VariableStruct.SumStartDay = xW(1,4)+1;
        end
    elseif strcmp(Type,'popupDates4W')
        if m_d(xW(1,3),2) == Date
            gen.VariableStruct.SumStartMonth = xW(1,3) + 1;
            gen.VariableStruct.SumStartDay = 1;
            if gen.VariableStruct.SumStartMonth == 13
                gen.VariableStruct.SumStartMonth = 1;
            end
        else
            gen.VariableStruct.SumStartMonth = xW(1,3);
            gen.VariableStruct.SumStartDay = Date + 1;
        end
    end
    days=(1:m_d(xW(1,1),2))';
    if get(handlesc(xW(2,2)),'Value')>length(days)
        set(handlesc(xW(2,2)),'Value',days(end))
    end
    set(handlesc(xW(2,2)),'String',days);
    clear days
    days=(1:m_d(xW(1,3),2))';
    if get(handlesc(xW(2,4)),'Value')>length(days)
        set(handlesc(xW(2,4)),'Value',days(end))
    end
    set(handlesc(xW(2,4)),'String',days);
end
end%Ends function popupDates

function gen = radiobuttonGridsellback(gen,hObject,eventdata,handles)
handles = guihandles;
if strcmp(get(hObject,'Tag'),'None')
    gen.VariableStruct.SellBackRate = 0;
    gen.VariableStruct.SellBackPerc = 0;
    gen.VariableStruct.MinImportThresh = 0;
    set(handles.textEdit11,'String','Minimum Import (kW)');
    set(handles.SellbackAsPerc,'Value',0);
    set(handles.FixedRate,'Value',0);
    set(handles.maxSellback,'String','0');
elseif strcmp(get(hObject,'Tag'),'SellbackAsPerc')
    gen.VariableStruct.SellBackRate = -1;
    gen.VariableStruct.SellBackPerc = str2double(get(handles.editTariffs,'String'));
    gen.VariableStruct.MinImportThresh = -inf;
    set(handles.textEdit11,'String','Maximum Export (kW)');
    set(handles.None,'Value',0);
    set(handles.FixedRate,'Value',0);
    set(handles.maxSellback,'String','inf');
elseif strcmp(get(hObject,'Tag'),'FixedRate')
    gen.VariableStruct.SellBackRate = str2double(get(handles.editSellbackRate,'String'));
    gen.VariableStruct.SellBackPerc = 0;
    gen.VariableStruct.MinImportThresh = -inf;
    set(handles.textEdit11,'String','Maximum Export (kW)');
    set(handles.SellbackAsPerc,'Value',0);
    set(handles.None,'Value',0);
    set(handles.maxSellback,'String','inf');
end
end%ends function radiobuttonGridsellback

function gen = editTariffs(gen,hObject,eventdata,handles)
handles = guihandles;
gen.VariableStruct.SellBackRate = -1;
gen.VariableStruct.SellBackPerc = str2double(get(handles.editTariffs,'String'));
end%Ends function editTariffs

function gen = editSellbackRate(gen,hObject,eventdata,handles)
handles = guihandles;
gen.VariableStruct.SellBackRate = str2double(get(handles.editSellbackRate,'String'));
gen.VariableStruct.SellBackPerc = 0;
end%ends function editSellbackRate

function gen = maxSellback(gen,hObject,eventdata,handles)
handles = guihandles;
if get(handles.None,'Value') == 0
    gen.VariableStruct.MinImportThresh = - str2double(get(handles.editMaxSellback,'String'));
else
    gen.VariableStruct.MinImportThresh = str2double(get(handles.editMaxSellback,'String'));
end
end%ends function maxSellback

function gen = radioSource(gen,hObject,eventdata,handles)
tags = {'NatGas';'Diesel';'Electricity';'Heat';};
for i =1:length(tags)
    if strcmp(tags{i},get(hObject,'String'))
        set(handles.(tags{i}),'Value',1)
    else
        set(handles.(tags{i}),'Value',0)
    end
end
if strcmp(get(hObject,'String'),'Natural Gas')
    gen.Source = 'NG';
else
    gen.Source = get(hObject,'String');
end
end%ends function radioSource

function gen = popupSolar(gen,hObject,eventdata,handles)
stateName = {'Alabama';'Alaska';'Arizona';'Arkansas';'California';'Colorado';'Connecticut';'Delaware';'Florida';'Georgia';'Hawaii';'Idaho';'Illinois';'Indiana';'Iowa';'Kansas';
         'Kentucky';'Louisiana';'Maine';'Maryland';'Massachusetts';'Michigan';'Minnesota';'Mississippi';'Missouri';'Montana';'Nebraska';'Nevada';'New Hampshire';'New Jersey';
         'New Mexico';'New York';'North Carolina';'North Dakota';'Ohio';'Oklahoma';'Oregon';'Pennsylvania';'Rhode Island';'South Carolina';'South Dakota';'Tennessee';'Texas';
         'Utah';'Vermont';'Virginia';'Washington';'West Virginia';'Wisconsin';'Wyoming';};
stateNum = get(hObject,'Value');
state = stateName(stateNum);
gen.VariableStruct.State = char(state);
end%ends function popupSolar

function gen = radiobuttonSType(gen,hObject,eventdata,handles)
tags = {'Flat Panel';'Concentrated';};
for i =1:length(tags)
    if strcmp(tags{i},get(hObject,'String'))
        set(handles.(tags{i}),'Value',1)
    else
        set(handles.(tags{i}),'Value',0)
    end
end
if strcmp(get(hObject,'String'),'Flat Panel')
    gen.VariableStruct.PVtype = 'flat';
elseif strcmp(get(hObject,'String'),'Concentrated')
    gen.VariableStruct.PVtype = 'concentrated';
end
end%ends function radiobuttonSType

function gen = radiobuttonSTracking(gen,hObject,eventdata,handles)
tags = {'Fixed';'SingleAxis';'DualAxis';};
for i =1:length(tags)
    if strcmp(tags{i},get(hObject,'String'))
        set(handles.(tags{i}),'Value',1)
    else
        set(handles.(tags{i}),'Value',0)
    end
end
if strcmp(get(hObject,'String'),'Fixed')
    gen.VariableStruct.Tracking = 'fixed';
elseif strcmp(get(hObject,'String'),'Single Axis')
    gen.VariableStruct.Tracking = '1axis';
elseif strcmp(get(hObject,'String'),'Dual Axis')
    gen.VariableStruct.Tracking = '2axis';
end
end%ends function radiobuttonSTracking

function gen = uitableDCAC(gen,hObject,eventdata,handles)
gen.VariableStruct.Data = get(hObject,'Data');
end%ends function uitableDCAC

function gen = uitableBat(gen,hObject,eventdata,handles)
gen.VariableStruct.VoltCurve = get(hObject,'Data');
end%ends function uitableBat

function gen = checkboxSeasonal(gen,hObject,eventdata,handles)
handles = guihandles;
str = get(hObject,'Tag');
switch str
    case 'sum_rate'
        set(handles.sum_rate,'Value',1);
        set(handles.win_rate,'Value',0);
        set(handles.compText7,'String',num2str(gen.VariableStruct.SumStartMonth));
        set(handles.compText8,'String',num2str(gen.VariableStruct.SumStartDay));
        set(handles.uitableEffCurve,'Data',gen.VariableStruct.SumRateTable);
        set(handles.compText1,'String',num2str(gen.VariableStruct.SumRates(1,1)));
        set(handles.compText2,'String',num2str(gen.VariableStruct.SumRates(1,2)));
        set(handles.compText3,'String',num2str(gen.VariableStruct.SumRates(2,1)));
        set(handles.compText4,'String',num2str(gen.VariableStruct.SumRates(2,2)));
        set(handles.compText5,'String',num2str(gen.VariableStruct.SumRates(3,1)));
        set(handles.compText6,'String',num2str(gen.VariableStruct.SumRates(3,2)));
    case 'win_rate'
        set(handles.sum_rate,'Value',0);
        set(handles.win_rate,'Value',1);
        set(handles.compText7,'String',num2str(gen.VariableStruct.WinStartMonth));
        set(handles.compText8,'String',num2str(gen.VariableStruct.WinStartDay));
        set(handles.uitableEffCurve,'Data',gen.VariableStruct.WinRateTable);
        set(handles.compText1,'String',num2str(gen.VariableStruct.WinRates(1,1)));
        set(handles.compText2,'String',num2str(gen.VariableStruct.WinRates(1,2)));
        set(handles.compText3,'String',num2str(gen.VariableStruct.WinRates(2,1)));
        set(handles.compText4,'String',num2str(gen.VariableStruct.WinRates(2,2)));
        set(handles.compText5,'String',num2str(gen.VariableStruct.WinRates(3,1)));
        set(handles.compText6,'String',num2str(gen.VariableStruct.WinRates(3,2)));
end
end%ends function checkboxSeasonal