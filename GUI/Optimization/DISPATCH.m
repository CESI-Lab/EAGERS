function varargout = DISPATCH(varargin)
% DISPATCH MATLAB code for DISPATCH.fig
%      DISPATCH, by itself, creates a new DISPATCH or raises the existing
%      singleton*.
%
%      H = DISPATCH returns the handle to a new DISPATCH or the handle to
%      the existing singleton*.
%
%      DISPATCH('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DISPATCH.M with the given input arguments.
%
%      DISPATCH('Property','Value',...) creates a new DISPATCH or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before DISPATCH_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to DISPATCH_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
 
% Edit the above text to modify the response to help DISPATCH
 
% Last Modified by GUIDE v2.5 28-Feb-2018 00:23:43
 
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @DISPATCH_OpeningFcn, ...
                   'gui_OutputFcn',  @DISPATCH_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end
 
if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT
 
 
% --- Executes just before DISPATCH is made visible.
function DISPATCH_OpeningFcn(hObject, eventdata, handles, varargin)
% Choose default command line output for EAGERS
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

global Plant GENINDEX mainFig TestData
GENINDEX = 1;
mainFig = gcf;
set(gcf,'Name','DISPATCH')
movegui(gcf,'center');

%% Set up Main Tabs
% Assumptions:
% 1. Tags of main tab static text boxes are of form, 'MainTab1',
% 'MainTab2', etc.
% 2. Tags of main tab panels are of form, 'uipanelMain1', 'uipanelMain2',
% etc.
TabText = {'Main Window';'Market Services';'Historian/Forecast';'Network View';'Settings'};
for i = 1:length(TabText)
    j = num2str(i);
    % panel management
    set(handles.(strcat('MainTab',j)),'Units','characters','String',TabText{i});
    set(handles.(strcat('uipanelMain',j)),'Units','characters','BorderType','none')
    if i ==1
        pan1pos = get(handles.uipanelMain1,'Position');
        pos = get(handles.MainTab1','Position');
        set(handles.MainTab1,'Position',[pos(1),pos(2),pos(3),pos(4)+.5])
    else
        pos2 = get(handles.(strcat('uipanelMain',j)),'Position');
        set(handles.(strcat('uipanelMain',j)),'Position',[pan1pos(1), pan1pos(2)+(pan1pos(4)-pos2(4)), pos2(3), pos2(4)])
        set(handles.(strcat('uipanelMain',j)),'Visible','off')
    end
end
set(handles.GenList,'Position',[1,16,40,27])
if ~isfield(Plant,'Building')
    Plant.Building = [];
end
[Plant.Generator,n_g] = check_ac_dc(Plant.Generator,Plant.Building,TestData);
n_b = length(Plant.Building);
list=cell(n_g+n_b,1);
for i=1:1:n_g
    list(i) = {Plant.Generator(i).Name};
end
for i=1:1:n_b
    list(n_g+i) = {Plant.Building(i).Name};
end
set(handles.uipanelMain1,'UserData',list)

network_names = fieldnames(Plant.Network);
network_names = network_names(~strcmp('name',network_names));
network_names = network_names(~strcmp('Equipment',network_names));
network_names = network_names(~strcmp('Location',network_names));
network_names = network_names(~strcmp('Buildings',network_names));
if isfield(Plant,'Building') && ~isempty(Plant.Building)
    network_names(end+1) = {'BuildingTemp'};
end
handles = PlotAxes_Setup(hObject, eventdata, handles,network_names,1);
handles = PlotAxes_Setup(hObject, eventdata, handles,network_names,3);

set(handles.popupmenuColor,'String',{'parula';'autumn';'cool';'spring';'summer';'winter';},'Value',1);

colormap(handles.ResultPlot1,'parula');
set(hObject,'UserData',colormap(handles.ResultPlot1));

handles = GenList_Make(handles);

set(handles.MarketPopup,'String',{'Market Prices';'Reserve Capacity';'Generation ($/kW)';'Demand Response ($/kW)'},'Value',1);
days = ceil(TestData.Timestamp(end) - TestData.Timestamp(1));
if days<7
    set(handles.sliderZoom,'Min',0,'Max',1,'Value',0,'SliderStep',[1,1]) %either single day or all data   
elseif days<31
    set(handles.sliderZoom,'Min',0,'Max',2,'Value',0,'SliderStep',[1/2,1/2]) %either single day, week or all data
elseif days<367
    set(handles.sliderZoom,'Min',0,'Max',3,'Value',0,'SliderStep',[1/3,1/3]) %either single day, week, month, or all data
else
    set(handles.sliderZoom,'Min',0,'Max',4,'Value',0,'SliderStep',[1/4,1/4]) %either single day, week, month, year, or all data
end
if days>20
    set(handles.sliderDate,'Min',0,'Max',days,'Value',0,'SliderStep',[1/days,10/days])
    set(handles.sliderStartDate,'Min',0,'Max',days,'Value',0,'SliderStep',[1/days,10/days])
else
    set(handles.sliderDate,'Min',0,'Max',days,'Value',0,'SliderStep',[1/days,1/days])
    set(handles.sliderStartDate,'Min',0,'Max',days,'Value',0,'SliderStep',[1/days,1/days])
end
insertMockups(handles)
switchTab(handles)

% --- Outputs from this function are returned to the command line.
function varargout = DISPATCH_OutputFcn(hObject, eventdata, handles) 
%ask user to save?

function switchTab(handles)
global GENINDEX 
tab = 1;
% Find out which tab is currently selected
while ~strcmp(get(handles.(strcat('uipanelMain',num2str(tab))),'Visible'),'on') %isempty(m) && isfield(handles,strcat('uipanelMain',num2str(i))) && 
    tab = tab+1;
end
if tab == 1
    set(handles.GenList,'Visible','on')
    if get(handles.Start,'Value') == 0
        slider_callback([])
    end
elseif tab == 2
    set(handles.GenList,'Visible','off')
    if get(handles.Start,'Value') == 0
        MarketServices(handles)
    end
elseif tab == 3
    GENINDEX = [];
    set(handles.GenList,'Visible','on')
    set(handles.Demands,'Visible','off')
    ForecastPlot
elseif tab == 4
    set(handles.GenList,'Visible','off')
elseif tab == 5
    set(handles.GenList,'Visible','on')
    loadSettingsTab(handles)
end

% Main tabs callback
function mainTab_Callback(hObject, eventdata, handles)
n = get(hObject,'Tag');
n = n(end);
i = 1;
% Find out which tab is currently selected
while ~strcmp(get(handles.(strcat('uipanelMain',num2str(i))),'Visible'),'on') %isempty(m) && isfield(handles,strcat('uipanelMain',num2str(i))) && 
    i = i+1;
end
m = num2str(i);

% CRUCIAL IN NEXT 3 STEPS: m, then n.
% Change color
bColor = get(handles.(strcat('MainTab',m)),'BackgroundColor');
set(handles.(strcat('MainTab',m)),'BackgroundColor',max(0,bColor-.1))
bColor = get(handles.(strcat('MainTab',n)),'BackgroundColor');
set(handles.(strcat('MainTab',n)),'BackgroundColor',min(1,bColor+.1))

% Change dimensions
pos = get(handles.(strcat('MainTab',m)),'Position');
set(handles.(strcat('MainTab',m)),'Position',[pos(1),pos(2),pos(3),pos(4)-.5])
pos = get(handles.(strcat('MainTab',n)),'Position');
set(handles.(strcat('MainTab',n)),'Position',[pos(1),pos(2),pos(3),pos(4)+.5])

%Change font color/weight
set(handles.(strcat('MainTab',m)),'FontWeight','normal','ForegroundColor',[0.501960784313726,0.501960784313726,0.501960784313726])
set(handles.(strcat('MainTab',n)),'FontWeight','bold','ForegroundColor',[0,0,0])

% Change visibility
set(handles.(strcat('uipanelMain',m)),'Visible','off')
set(handles.(strcat('uipanelMain',n)),'Visible','on')
switchTab(handles)

function insertMockups(handles)
global Model_dir
handles.NetworkDisplay = axes('Units','normalized',...
        'Position', [0,.1,1,.8],...
        'Tag', 'NetworkDisplay',...
        'Parent', handles.uipanelMain4,...
        'Visible','on');
[x,map] = imread(fullfile(Model_dir,'GUI','Graphics','NetworkDisplayMockup.png'));
image(x,'Parent',handles.NetworkDisplay);
set(handles.NetworkDisplay,'xtick',[],'xticklabel',[],'ytick',[],'yticklabel',[],'box','off')

function loadSettingsTab(handles)
global Plant
%Set handles for Settings
set(handles.editCommandOnOff,'string','--')
set(handles.editCommandSet,'string','--')
set(handles.editMeasureOnOff,'string','--')
set(handles.editMeasureInput,'string','--')
set(handles.editMeasurePrimary,'string','--')
set(handles.editMeasureSecondary,'string','--')

set(handles.Interval,'string',Plant.optimoptions.Interval);
set(handles.Horizon, 'string', Plant.optimoptions.Horizon);
set(handles.Resolution, 'string', Plant.optimoptions.Resolution);
set(handles.SimBuilding,'value',Plant.optimoptions.EnergyPlus);

set(handles.uipanelSolverMethod, 'Visible', 'on');
set(handles.changingtimesteps, 'Visible', 'on');
set(handles.OtherOptions, 'Visible', 'on');
set(handles.constant, 'value', strcmp(Plant.optimoptions.tspacing,'constant'));
set(handles.linear, 'value', strcmp(Plant.optimoptions.tspacing, 'linear'));
set(handles.logarithm, 'value', strcmp(Plant.optimoptions.tspacing, 'logarithm'));
set(handles.manual, 'value', strcmp(Plant.optimoptions.tspacing, 'manual'));
set(handles.scaletime, 'string', Plant.optimoptions.scaletime);
set(handles.excessHeat, 'value', Plant.optimoptions.excessHeat);
set(handles.excessCool, 'value', Plant.optimoptions.excessCool);
if strcmp(Plant.optimoptions.solver,'NREL')
    set(handles.nrel_gurobi, 'value', true);
    set(handles.changingtimesteps, 'Visible', 'off');
    set(handles.OtherOptions, 'Visible', 'off');
elseif strcmp(Plant.optimoptions.solver,'ANN')
    set(handles.ANNradiobutton, 'value', true);
else
    set(handles.NoMixedInteger, 'value', ~Plant.optimoptions.MixedInteger);
    set(handles.MixedInteger, 'value', Plant.optimoptions.MixedInteger);
end
set(handles.SpinReserve, 'value', Plant.optimoptions.SpinReserve);
set(handles.SpinReservePerc, 'string', Plant.optimoptions.SpinReservePerc);
set(handles.Topt, 'string', Plant.optimoptions.Topt);
set(handles.Tmpc, 'string', Plant.optimoptions.Tmpc);

n_g = length(Plant.Generator);
str = {};
stor = [];
for i = 1:1:n_g
    if ismember(Plant.Generator(i).Type,{'Electric Storage';'Thermal Storage'; 'Hydro Storage';})
        str{end+1} = Plant.Generator(i).Name;
        stor(end+1) = i;
        if ~isfield(Plant.Generator(i).VariableStruct,'Buffer')
            Plant.Generator(i).VariableStruct.Buffer = 20;
        end
    end
end
if ~isempty(str)
    set(handles.StorageBuff, 'Visible', 'on','UserData',stor);
    set(handles.editBuffer, 'Visible', 'on');
    set(handles.textBuffer, 'Visible', 'on');
    set(handles.StorageBuff,'string',str,'value',1);
    set(handles.editBuffer, 'string', Plant.Generator(stor(1)).VariableStruct.Buffer);
else
    set(handles.StorageBuff, 'Visible', 'off');
    set(handles.editBuffer, 'Visible', 'off');
    set(handles.textBuffer, 'Visible', 'off');
end

if isfield(Plant,'Building') && ~isempty(Plant.Building)
%     Plant.optimoptions.forecast ='Building';
%     set(handles.ForecastMethod,'Visible','off');
    set(handles.Building, 'Visible','on');
else
    set(handles.Building, 'Visible','off');
end
set(handles.ForecastMethod,'Visible','on');
set(handles.ARMA, 'value', false);
set(handles.ARIMA, 'value', false);
set(handles.NeuralNet, 'value', false);
set(handles.Surface, 'value', false);
set(handles.Perfect, 'value', false);
set(handles.Building, 'value',false);
set(handles.(Plant.optimoptions.forecast), 'value', true);

 
% --- Executes on button press in Switch.
function Switch_Callback(hObject, eventdata, handles)
global mainFig
mainFig = [];
Stop_Callback([],[],handles);
%%send user back to EPT, pass along the plant generators 
close
MainScreen1

% --- Creates plots on the Main Window and Forecasting Tabs
function handles = PlotAxes_Setup(hObject, eventdata, handles,networkNames,tab)
nPlot = length(networkNames);
if nPlot ==1 %1 very large plot
    Pos = [52 6 150 30];
    vis = {'on'};
elseif nPlot<=4 %1 large plot and n-1 smaller plots
    Pos = [52 10 80 26;];
    vis = {'on'};
    for j = 2:nPlot
        Pos(j,:) = [153 3+(j-2)*40/(nPlot-1) 50 min(20,(40/(nPlot-1)-6))];
        vis(end+1) = {'on'};
    end
else%plot three at a time
    Pos = [52 10 80 26;153 3 50 14;153 23 50 14];
    vis = {'on';'on';'on';};
    for j = 4:1:nPlot
        Pos(end+1,:) = [153 3+20*rem(j,2) 50 14];
        vis(end+1) = {'off'};
    end
end
if nPlot>4
    set(handles.NextPlot,'Visible','on','UserData',1)
    set(handles.PrevPlot,'Visible','on','UserData',1)
else
    set(handles.NextPlot,'Visible','off')
    set(handles.PrevPlot,'Visible','off')
end
if tab==1
    name = 'Result';
else
    name = 'Forecast';
end
quote='''';
for i = 1:1:nPlot
    %primary axes
    handles.(strcat(name,'Plot',num2str(i))) = axes('Units','characters',...
        'Position', Pos(i,:),'NextPlot','add',...
        'Tag', strcat(name,'Plot',num2str(i)),...
        'Parent', handles.(strcat('uipanelMain',num2str(tab))),...
        'Visible',vis{i});
    %secondary axes (y-label on right)
    handles.(strcat(name,'Plot',num2str(i),'b')) = axes('Units','characters',...
        'Position', Pos(i,:),'NextPlot','add','color','none',...
        'Tag', strcat(name,'Plot',num2str(i),'b'),...
        'Parent', handles.(strcat('uipanelMain',num2str(tab))),...
        'YAxisLocation','right',...
        'Visible',vis{i});
    callback = strcat('@(hObject,eventdata)DISPATCH(',quote,'PlotSwitch',quote,',hObject)');
    handles.(strcat(name,'Name',num2str(i))) = uicontrol('Style', 'pushbutton', 'String', networkNames{i},...
        'Units','characters',...
        'Position', [Pos(i,1)+Pos(i,3)/10,Pos(i,2)+Pos(i,4)+.25,35,2],...
        'BackGroundColor',[0.94 0.94 0.94],...
        'Tag', strcat(name,'Name',num2str(i)),...
        'FontSize',16,...
        'FontUnits','points',...
        'FontWeight','bold',...
        'Parent', handles.(strcat('uipanelMain',num2str(tab))),...
        'Callback',eval(callback),...
        'Visible',vis{i},...
        'UserData',i);
    if tab==1
        callback = strcat('@(hObject,eventdata)DISPATCH(',quote,'show_different_storage_Callback',quote,',hObject)');
        handles.(strcat(name,'StorageSelect',num2str(i))) = uicontrol('Style', 'popupmenu', 'String', networkNames{i},...
            'Units','characters',...
            'Position', [Pos(i,1)+40,Pos(i,2)+Pos(i,4)+.25,30,2],...
            'Tag', strcat(name,'StorageSelect',num2str(i)),...
            'FontSize',11,...
            'FontUnits','points',...
            'FontWeight','normal',...
            'Parent', handles.(strcat('uipanelMain',num2str(tab))),...
            'Callback',eval(callback),...
            'Visible','off');
    end
end

function PlotSwitch(hObject)
%if you pick the main plot and it is network1, swap with plot 2, 
% otherwise swap the small plot picked and main plot
global Plant
handles = guihandles;
networkNames = fieldnames(Plant.subNet);
nPlot = length(networkNames);
net = get(hObject,'String');
plotI = get(hObject,'UserData');
netI = find(strcmp(net,networkNames));
if strcmp(get(handles.uipanelMain1,'Visible'),'on')
    screen = 'Result';
elseif strcmp(get(handles.uipanelMain3,'Visible'),'on')
    screen = 'Forecast';
end
for i = 1:1:nPlot %put names back to default pos
    set(handles.(strcat(screen,'Name',num2str(i))),'String',networkNames{i});
end
    
if netI == plotI %graph is currently in its default position (move to primary)
    if netI == 1 && nPlot>1%already is primary swap with 2
        set(handles.(strcat(screen,'Name1')),'String',networkNames{2});
        set(handles.(strcat(screen,'Name2')),'String',networkNames{1});
    else
        set(handles.(strcat(screen,'Name1')),'String',networkNames{netI});
        set(handles.(strcat(screen,'Name',num2str(netI))),'String',networkNames{1});
    end
end
if strcmp(screen,'Forecast')
    ForecastPlot
elseif strcmp(screen,'Result')
    slider_callback([])
end


% --- Executes on button press in PrevPlot.
function PrevPlot_Callback(hObject, eventdata, handles)
global Plant
handles = guihandles;
networkNames = fieldnames(Plant.subNet);
nPlot = length(networkNames);
if strcmp(get(handles.uipanelMain1,'Visible'),'on')
    screen = 'Result';
elseif strcmp(get(handles.uipanelMain3,'Visible'),'on')
    screen = 'Forecast';
end
net = get(handles.(strcat(screen,'Name2')),'String');
netI = find(strcmp(net,networkNames));
if netI == 1
    netI = nPlot;
else
    netI = netI - 1;
end
set(handles.(strcat(screen,'Name2')),'String',networkNames{netI});

net = get(handles.(strcat(screen,'Name3')),'String');
netI = find(strcmp(net,networkNames));
if netI == 1
    netI = nPlot;
else
    netI = netI - 1;
end
set(handles.(strcat(screen,'Name3')),'String',networkNames{netI});
if strcmp(screen,'Forecast')
    ForecastPlot
elseif strcmp(screen,'Result')
    slider_callback([])
end

% --- Executes on button press in NextPlot.
function NextPlot_Callback(hObject, eventdata, handles)
global Plant
handles = guihandles;
networkNames = fieldnames(Plant.subNet);
nPlot = length(networkNames);
if strcmp(get(handles.uipanelMain1,'Visible'),'on')
    screen = 'Result';
elseif strcmp(get(handles.uipanelMain3,'Visible'),'on')
    screen = 'Forecast';
end
net = get(handles.(strcat(screen,'Name2')),'String');
netI = find(strcmp(net,networkNames));
if netI == nPlot
    netI = 1;
else
    netI = netI + 1;
end
set(handles.(strcat(screen,'Name2')),'String',networkNames{netI});

net = get(handles.(strcat(screen,'Name3')),'String');
netI = find(strcmp(net,networkNames));
if netI == nPlot
    netI = 1;
else
    netI = netI + 1;
end
set(handles.(strcat(screen,'Name3')),'String',networkNames{netI});
if strcmp(screen,'Forecast')
    ForecastPlot
elseif strcmp(screen,'Result')
    slider_callback([])
end

function handles = GenList_Make(handles)
global Plant Model_dir
quote='''';
list = get(handles.uipanelMain1,'UserData');
if length(list)>12
    set(handles.NextGen,'Visible','on','UserData',1)
    set(handles.PrevGen,'UserData',1)
end
%Makes buttons for GUI && status buttons next to corresponding generator
nG = length(Plant.Generator);
if isfield(Plant,'Building') && ~isempty(Plant.Building)
    nB = length(Plant.Building);
else
    nB = 0;
end
colorVec = get(gcf,'UserData');
colorsPlot = interp1(linspace(0,1,length(colorVec)),colorVec,linspace(0,1,nG+nB));
for i=1:1:nG+nB
    curtab = floor((i-1)/12)+1;
    prev = 12*(curtab-1);
    if curtab==1
        vis = 'on';
    else
        vis = 'off';
    end
    pos = 27 - 2*(i-prev);
    if i<=nG
        num = num2str(i);
        callback = strcat('@(hObject,eventdata)DISPATCH(',quote,'Gen_Callback',quote,',hObject,eventdata,guidata(hObject))');
        handles.(strcat('Generator_',num)) = uicontrol('Style', 'pushbutton', 'String', list{i},...
        'Units','characters','Position', [1 pos 25 1.8],...
        'Tag', strcat('Generator_',num),'FontSize', 10,...
        'Parent', handles.GenList,...
        'Callback',eval(callback),'Visible',vis,'UserData',i,'BackgroundColor',colorsPlot(i,:));
    else
        num = num2str(i-nG);
        callback = strcat('@(hObject,eventdata)DISPATCH(',quote,'Build_Callback',quote,',hObject,eventdata,guidata(hObject))');
        handles.(strcat('Building_',num)) = uicontrol('Style', 'pushbutton', 'String', list{i},...
        'Units','characters','Position', [1 pos 25 1.8],...
        'Tag', strcat('Building_',num),'FontSize', 10,...
        'Parent', handles.GenList,...
        'Callback',eval(callback),'Visible',vis,'UserData',i-nG,'BackgroundColor',colorsPlot(i,:));
    end
    if i<=nG %Only make Status buttons on Main Window
        pos = 27.5 - 2*(i-prev);
        callback = strcat('@(hObject,eventdata)DISPATCH(',quote,'Status_Callback',quote,',hObject,eventdata,guidata(hObject))');
        if Plant.Generator(i).Enabled
            [x,map] = imread(fullfile(Model_dir,'GUI','Graphics','green.png'));
        else
            [x,map] = imread(fullfile(Model_dir,'GUI','Graphics','red.png'));
        end
        if license('test','Image Processing Toolbox')
            set(handles.Switch,'Units','pixels');
            posP = get(handles.Switch,'Position');
            set(handles.Switch,'Units','characters');
            posC = get(handles.Switch,'Position');
            x = imresize(x,[3*posP(3)/posC(3) posP(4)/posC(4)]);
        end
        if Plant.Generator(i).Enabled
            enableGen  = 'bold';
        else
            enableGen  = 'normal';
        end
        handles.(strcat('GeneratorStat_',num)) = uicontrol('Style', 'pushbutton', 'String', '',...
        'Units','characters',...
        'Position', [27 pos 3 1],...
        'Tag', strcat('GeneratorStat_',num),...
        'cdata', x,...
        'FontWeight',enableGen,...
        'Parent', handles.GenList,...
        'Callback',eval(callback),...
        'Visible',vis,...
        'UserData',i);
    end
end
% --- Executes on button press in PrevGen.
function PrevGen_Callback(hObject, eventdata, handles)
handles = guihandles;
list = get(handles.uipanelMain1,'UserData');
gen = length(list);
page = get(handles.PrevGen,'UserData');%current page of the list
if (page-1)<2%if the new page is the 1st
    set(handles.PrevGen,'Visible','off','UserData',page-1)
else
    set(handles.PrevGen,'Visible','on','UserData',page-1);
end
set(handles.NextGen,'Visible','on','UserData',page-1)
for i = 1:1:12
    if 12*(page-1)+i<=gen
         j = num2str(12*(page-1)+i);
         set(handles.(strcat('Generator_',j)),'Visible','off');
         set(handles.(strcat('GeneratorStat_',j)),'Visible','off');
    end
end
for i = 1:1:12
    j = num2str(12*(page-2)+i);
     set(handles.(strcat('Generator_',j)),'Visible','on');
     set(handles.(strcat('GeneratorStat_',j)),'Visible','on');
end


% --- Executes on button press in NextGen.
function NextGen_Callback(hObject, eventdata, handles)
handles = guihandles;
list = get(handles.uipanelMain1,'UserData');
gen = length(list);
page = get(handles.PrevGen,'UserData');%current page of the list
if page==ceil(gen/12)-1
    set(handles.NextGen,'Visible','off','UserData',page+1)
else
    set(handles.NextGen,'Visible','on','UserData',page+1);
end
set(handles.PrevGen,'Visible','on','UserData',page+1)
for i = 1:1:12
     j = num2str(12*(page-1)+i);
     set(handles.(strcat('Generator_',j)),'Visible','off');
     set(handles.(strcat('GeneratorStat_',j)),'Visible','off');
end
for i = 1:1:12
    j = num2str(12*(page)+i);
    if 12*(page)+i<=gen
         set(handles.(strcat('Generator_',j)),'Visible','on');
         set(handles.(strcat('GeneratorStat_',j)),'Visible','on');
    end
end

%When Component Buttons on the main tab are clicked
function Gen_Callback(hObject, eventdata, handles)
global Plant GENINDEX BUILDINDEX
tab = 1;
% Find out which tab is currently selected
while ~strcmp(get(handles.(strcat('uipanelMain',num2str(tab))),'Visible'),'on') %isempty(m) && isfield(handles,strcat('uipanelMain',num2str(i))) && 
    tab = tab+1;
end
if tab == 1
    GENINDEX = get(hObject,'UserData');
    size = num2str(Plant.Generator(GENINDEX).Size);
    set(handles.SelGen,'Title',Plant.Generator(GENINDEX).Name,'UserData',GENINDEX)
    if strcmp(Plant.Generator(GENINDEX).Type,'Utility')
        set(handles.GenSpec1,'String','Inf')
    else
        set(handles.GenSpec1,'String',size)
    end
    if Plant.Generator(GENINDEX).Enabled == 1
        set(handles.GenEnable,'Value',1)
        set(handles.GenDisable,'Value',0)
    elseif Plant.Generator(GENINDEX).Enabled == 0
        set(handles.GenEnable,'Value',0)
        set(handles.GenDisable,'Value',1)
    end
elseif tab == 3
    %plot the forecastmethod and actual dispatch of a single generator
    BUILDINDEX = [];
    GENINDEX = get(hObject,'UserData');
    set(handles.Demands,'Visible','on')
    ForecastPlot
elseif tab == 5
    i = get(hObject,'UserData');
    set(handles.CurrentComm,'String',['Ports for: ',Plant.Generator(i).Name])
    if isfield(Plant.Generator(i).VariableStruct,'Comm')
        Comm = Plant.Generator(i).VariableStruct.Comm;
    else
        Comm.OnOff = 0;
        Comm.Set = 0;
    end
    if isfield(Plant.Generator(i).VariableStruct,'Measure')
        Measure = Plant.Generator(i).VariableStruct.Measure;
    else
        Measure.OnOff = 0;
        Measure.Input = 0;
        Measure.Electric = 0;
        Measure.Thermal = 0;
    end
    set(handles.editCommandOnOff,'string',num2str(Comm.OnOff))
    set(handles.editCommandSet,'string',num2str(Comm.Set))
    set(handles.editMeasureOnOff,'string',num2str(Measure.OnOff))
    set(handles.editMeasureInput,'string',num2str(Measure.Input))
    set(handles.editMeasurePrimary,'string',num2str(Measure.Electric))
    set(handles.editMeasureSecondary,'string',num2str(Measure.Thermal))
end

%When Building Buttons on the main tab are clicked
function Build_Callback(hObject, eventdata, handles)
global Plant GENINDEX BUILDINDEX
tab = 1;
% Find out which tab is currently selected
while ~strcmp(get(handles.(strcat('uipanelMain',num2str(tab))),'Visible'),'on') %isempty(m) && isfield(handles,strcat('uipanelMain',num2str(i))) && 
    tab = tab+1;
end
if tab == 1
    b = get(hObject,'UserData');
    build = Plant.Building(b);
elseif tab ==3
    GENINDEX = [];
    BUILDINDEX = get(hObject,'UserData');
    set(handles.Demands,'Visible','on')
    ForecastPlot
elseif tab ==5
    i = get(hObject,'UserData');
    set(handles.CurrentComm,'String',['Ports for: ',Plant.Building(i).Name])
    if isfield(Plant.Building(i).VariableStruct,'Measure')
        Measure = Plant.Building(i).VariableStruct.Measure;
    else
        Measure.Temperature = 0;
        Measure.Humidity = 0;
    end
    set(handles.editCommandOnOff,'visible','off')
    set(handles.editCommandSet,'visible','off')
    set(handles.editMeasureOnOff,'visible','off')
    set(handles.editMeasureInput,'visible','off')
    set(handles.editMeasurePrimary,'string',num2str(Measure.Temperature))
    set(handles.editMeasureSecondary,'string',num2str(Measure.Humidity))
end

%When Status colors are clicked
function Status_Callback(hObject, eventdata,handles)
global Plant Model_dir
i = get(hObject,'UserData');
size = num2str(Plant.Generator(i).Size);
if ~isempty(strfind(Plant.Generator(i).Name,'Utility'))
    set(handles.GenSpec1,'String','Inf')
else
    set(handles.GenSpec1,'String',size)
end
set(handles.SelGen,'Title',Plant.Generator(i).Name,'UserData',i)
if strcmp(get(hObject,'FontWeight'),'bold')
    Plant.Generator(i).Enabled = 0;
    set(hObject,'FontWeight','normal')
    [x,map] = imread(fullfile(Model_dir,'GUI','Graphics','red.png'));
else
    Plant.Generator(i).Enabled = 1;
    set(hObject,'FontWeight','bold')
    [x,map] = imread(fullfile(Model_dir,'GUI','Graphics','green.png'));
end
set(handles.GenEnable,'value',Plant.Generator(i).Enabled)
set(handles.GenDisable,'Value',~Plant.Generator(i).Enabled)
if license('test','Image Processing Toolbox')
    set(handles.Switch,'Units','pixels');
    pos1 = get(handles.Switch,'Position');
    set(handles.Switch,'Units','characters');
    pos2 = get(handles.Switch,'Position');
    x = imresize(x,[3*pos1(3)/pos2(3) pos1(4)/pos2(4)]);
end
set(hObject,'cdata',x)


function sliderStartDate_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes on button press in Start.
function Start_Callback(hObject, eventdata, handles)
global Plant TestData  %Virtual RealTime 
%Virtual: Running a simulation only, set to zero when the end of the test data set is reached
% if get(handles.VirtualMode,'Value') ==1 
%     Virtual = 1;
%     RealTime = 0;
%     Plant.optimoptions.mode = 'virtual';
% elseif get(handles.ObserverMode,'Value') == 1
%     Virtual = 1;
%     RealTime = 0;
%     Plant.optimoptions.mode = 'observer';
% elseif get(handles.ControllerMode,'Value') == 1
%     Virtual = 0;
%     RealTime = 1;
%     Plant.optimoptions.mode = 'controller';
% end

handles = guihandles;
date = TestData.Timestamp(1) + get(handles.sliderStartDate,'Value');

set(handles.Start,'Value',1);%reset start button
set(handles.Stop,'Value',0);%reset stop button
if Plant.optimoptions.EnergyPlus
    %% Start EnergyPlus cosimulation
    dir = strrep(which('run_simulation.m'),fullfile('run_simulation.m'),'');
    cd(fullfile(dir,'NRELmethod'))
    installMlep
    ep = mlepProcess;
    %% need to automate the building and weather file selection
    ep.arguments = {'B_2ndCoSim_v2', Plant.Weather.filename};
    ep.acceptTimeout = 6000;
    [status, msg] = ep.start;  
    if status ~= 0
        error('Could not start EnergyPlus: %s.', msg);
    end
    [status, msg] = ep.acceptSocket;
    if status ~= 0
        error('Could not connect to EnergyPlus: %s.', msg);
    end
else
    ep = [];
end
if ~isfield(Plant,'Building')
    Plant.Building = [];
end
if strcmp(Plant.optimoptions.solver,'NREL')
    [Plant.Generator,Plant.Building,Plant.Dispatch,Plant.Predicted] = ...
        run_simulation_NREL(Plant.Generator,Plant.Building,TestData,Plant.optimoptions,date,handles,ep); %send to a different function that links to energy plus co-simulation
else
    if ~isfield(Plant,'fluid_loop') 
        Plant.fluid_loop = [];
    end
    data = [];
    if isfield(Plant,'Data') 
        data = Plant.Data;
    end
    weather = [];
    if isfield(Plant,'Data') 
        weather = Plant.Data;
    end
    TestData = update_test_data(TestData,data,weather,Plant.optimoptions);
    num_steps = 0;
    s_i = 1;
    if isfield(Plant,'Dispatch') && isfield(Plant.Dispatch,'Timestamp')
        if any(Plant.Dispatch.Timestamp)
            s_i = 1+floor((date-Plant.Dispatch.Timestamp(1))*24/Plant.optimoptions.Resolution);
            if s_i>1 && abs(Plant.Dispatch.Timestamp(s_i) - date)<1e-7 %resuming after sucessful hitting stop button
                num_steps = length(Plant.Dispatch.Timestamp); %number of simulation steps
            end
        end
    end
    if num_steps == 0
        [Plant.Generator,Plant.Building,Plant.fluid_loop,Plant.subNet,Plant.OpMatA,Plant.OpMatB,Plant.OneStep,~,num_steps,Plant.Dispatch,Plant.Predicted,Plant.Design,flag] = ...
            initialize_optimization(Plant.Generator,Plant.Building,Plant.fluid_loop,Plant.Network,Plant.optimoptions,TestData,date);
    end
    [Plant.Generator,Plant.Building,Plant.fluid_loop,Plant.Design,Plant.Dispatch,Plant.Predicted] = ...
        run_simulation(date,num_steps,s_i,handles,TestData,Plant.Generator,Plant.Building,Plant.fluid_loop,Plant.optimoptions,Plant.subNet,Plant.OpMatA,Plant.OpMatB,Plant.OneStep,Plant.Design,Plant.Dispatch,Plant.Predicted,ep);
    if isfield(Plant.Dispatch,'OutFlow')
        TestData.Hydro.OutFlow = Plant.Dispatch.OutFlow;
    end
end
Stop_Callback(hObject, eventdata, handles)

% --- Executes on button press in Stop.
function Stop_Callback(hObject, eventdata, handles)
%stops dispatch
global DispatchWaitbar %Virtual RealTime
set(handles.Stop,'Value',1);%reset stop button
set(handles.Start,'Value',0);%reset start button
% Virtual = 0;
% if get(handles.ControllerMode,'Value') == 1
%     if RealTime
%         closePorts;
%     end
%     RealTime = 0;%end condition for real simulation
%     Virtual = 0;%end condition for virtual simulation
% 
%     T1 = timerfind('Name', 'dispTimer') ;
%     T2 = timerfind('Name', 'optTimer') ;
%     T3 = timerfind('Name', 'mpcTimer') ;
%     T4 = timerfind('Name', 'fanTimer') ;
%     Timers = [T1,T2,T3,T4];
%     for i = 1:1:length(Timers)
%         stop(Timers(i));
%         delete(Timers(i))
%     end
% end
close(DispatchWaitbar)
DispatchWaitbar=[];

% --- Executes on button press in GenEnable.
function GenEnable_Callback(hObject, eventdata, handles)
global Plant Model_dir
i = get(handles.SelGen,'UserData');
Plant.Generator(i).Enabled = 1;
[x,map] = imread(fullfile(Model_dir,'GUI','Graphics','green.png'));
if license('test','Image Processing Toolbox')
    set(handles.Switch,'Units','pixels');
    pos1 = get(handles.Switch,'Position');
    set(handles.Switch,'Units','characters');
    pos2 = get(handles.Switch,'Position');
    x = imresize(x,[3*pos1(3)/pos2(3) pos1(4)/pos2(4)]);
end
num = 2*length(Plant.Generator) - (2*i - 1);
set(handles.uipanelMain1.Children(num),'FontWeight','normal','cdata',x)

% --- Executes on button press in GenDisable.
function GenDisable_Callback(hObject, eventdata, handles)
global Plant Model_dir
i = get(handles.SelGen,'UserData');
Plant.Generator(i).Enabled = 0;
[x,map] = imread(fullfile(Model_dir,'GUI','Graphics','red.png'));
if license('test','Image Processing Toolbox')
    set(handles.Switch,'Units','pixels');
    pos1 = get(handles.Switch,'Position');
    set(handles.Switch,'Units','characters');
    pos2 = get(handles.Switch,'Position');
    x = imresize(x,[3*pos1(3)/pos2(3) pos1(4)/pos2(4)]);
end
num = 2*length(Plant.Generator) - (2*i - 1);
set(handles.uipanelMain1.Children(num),'FontWeight','normal','cdata',x)
        
%% for manual control?
function GenStatus1_Callback(hObject, eventdata, handles)
function GenStatus1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function GenStatus2_Callback(hObject, eventdata, handles)
function GenStatus2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function GenSpec1_Callback(hObject, eventdata, handles)
function GenSpec1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function GenSpec2_Callback(hObject, eventdata, handles)
function GenSpec2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in LineGraph.
function LineGraph_Callback(hObject, eventdata, handles)
if get(handles.StackedGraph,'Value')==1
    a = get(handles.StackedGraph,'BackgroundColor');
    b = get(handles.LineGraph,'BackgroundColor');
    c = get(handles.StackedGraph,'ForegroundColor');
    d = get(handles.LineGraph,'ForegroundColor');
    set(handles.LineGraph,'Value',1,'BackgroundColor',a,'ForegroundColor',c);
    set(handles.StackedGraph,'Value',0,'BackgroundColor',b,'ForegroundColor',d);
else
    set(handles.LineGraph,'Value',1); %was already pressed
end
slider_callback([])


% --- Executes on button press in StackedGraph.
function StackedGraph_Callback(hObject, eventdata, handles)
if get(handles.LineGraph,'Value')==1
    a = get(handles.StackedGraph,'BackgroundColor');
    b = get(handles.LineGraph,'BackgroundColor');
    c = get(handles.StackedGraph,'ForegroundColor');
    d = get(handles.LineGraph,'ForegroundColor');
    set(handles.LineGraph,'Value',0,'BackgroundColor',a,'ForegroundColor',c);
    set(handles.StackedGraph,'Value',1,'BackgroundColor',b,'ForegroundColor',d);
else
    set(handles.StackedGraph,'Value',1); %was already pressed
end
slider_callback([])

function show_different_storage_Callback(hObject)
slider_callback([])

% --- Executes on button press in AutoControl.
function AutoControl_Callback(hObject, eventdata, handles)
if get(handles.AutoControl,'Value')==1
    a = get(handles.ManualControl,'BackgroundColor');
    b = get(handles.AutoControl,'BackgroundColor');
    c = get(handles.ManualControl,'ForegroundColor');
    d = get(handles.AutoControl,'ForegroundColor');
    set(handles.AutoControl,'Value',1,'BackgroundColor',a,'ForegroundColor',c);
    set(handles.ManualControl,'Value',0,'BackgroundColor',b,'ForegroundColor',d);
else
    set(handles.AutoControl,'Value',1); %was already pressed
end

% --- Executes on button press in ManualControl.
function ManualControl_Callback(hObject, eventdata, handles)
if get(handles.ManualControl,'Value')==1
    a = get(handles.ManualControl,'BackgroundColor');
    b = get(handles.AutoControl,'BackgroundColor');
    c = get(handles.ManualControl,'ForegroundColor');
    d = get(handles.AutoControl,'ForegroundColor');
    set(handles.AutoControl,'Value',0,'BackgroundColor',a,'ForegroundColor',c);
    set(handles.ManualControl,'Value',1,'BackgroundColor',b,'ForegroundColor',d);
else
    set(handles.ManualControl,'Value',1); %was already pressed
end

% --- Executes on button press in VirtualMode.
function VirtualMode_Callback(hObject, eventdata, handles)
if get(handles.VirtualMode,'Value')==1
    a = get(handles.VirtualMode,'BackgroundColor');
    c = get(handles.VirtualMode,'ForegroundColor');
    if get(handles.ObserverMode,'Value')==1
        b = get(handles.ObserverMode,'BackgroundColor');
        d = get(handles.ObserverMode,'ForegroundColor');
        set(handles.ObserverMode,'Value',0,'BackgroundColor',a,'ForegroundColor',c);
    elseif get(handles.ControllerMode,'Value')==1
        b = get(handles.ControllerMode,'BackgroundColor');
        d = get(handles.ControllerMode,'ForegroundColor');
        set(handles.ControllerMode,'Value',0,'BackgroundColor',a,'ForegroundColor',c);
    end
    set(handles.VirtualMode,'Value',1,'BackgroundColor',b,'ForegroundColor',d);
else
    set(handles.VirtualMode,'Value',1); %was already pressed
end

% --- Executes on button press in ObserverMode.
function ObserverMode_Callback(hObject, eventdata, handles)
if get(handles.ObserverMode,'Value')==1
    a = get(handles.ObserverMode,'BackgroundColor');
    c = get(handles.ObserverMode,'ForegroundColor');
    if get(handles.VirtualMode,'Value')==1
        b = get(handles.VirtualMode,'BackgroundColor');
        d = get(handles.VirtualMode,'ForegroundColor');
        set(handles.VirtualMode,'Value',0,'BackgroundColor',a,'ForegroundColor',c);
    elseif get(handles.ControllerMode,'Value')==1
        b = get(handles.ControllerMode,'BackgroundColor');
        d = get(handles.ControllerMode,'ForegroundColor');
        set(handles.ControllerMode,'Value',0,'BackgroundColor',a,'ForegroundColor',c);
    end
    set(handles.ObserverMode,'Value',1,'BackgroundColor',b,'ForegroundColor',d);
else
    set(handles.ObserverMode,'Value',1); %was already pressed
end

% --- Executes on button press in ControllerMode.
function ControllerMode_Callback(hObject, eventdata, handles)
if get(handles.ControllerMode,'Value')==1
    a = get(handles.ControllerMode,'BackgroundColor');
    c = get(handles.ControllerMode,'ForegroundColor');
    if get(handles.ObserverMode,'Value')==1
        b = get(handles.ObserverMode,'BackgroundColor');
        d = get(handles.ObserverMode,'ForegroundColor');
        set(handles.ObserverMode,'Value',0,'BackgroundColor',a,'ForegroundColor',c);
    else
        b = get(handles.VirtualMode,'BackgroundColor');
        d = get(handles.VirtualMode,'ForegroundColor');
        set(handles.VirtualMode,'Value',0,'BackgroundColor',a,'ForegroundColor',c);
    end
    set(handles.ControllerMode,'Value',1,'BackgroundColor',b,'ForegroundColor',d);
else
    set(handles.ControllerMode,'Value',1); %was already pressed
end

%% Historian/ForecastMethod Tab
function Demands_Callback(hObject, eventdata, handles)
%plot the forecastmethod and actual demands
global GENINDEX BUILDINDEX
GENINDEX = [];
BUILDINDEX = [];
set(handles.Demands,'Visible','off')
ForecastPlot

function sliderDate_Callback(hObject, eventdata, handles)
ForecastPlot
function sliderDate_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
function sliderZoom_Callback(hObject, eventdata, handles)
ForecastPlot
function sliderZoom_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function ForecastPlot
global Plant GENINDEX TestData
handles = guihandles;
%find the current date
if ~isempty(GENINDEX)
    slider_callback(GENINDEX)
else
    date_1 = TestData.Timestamp(1) + get(handles.sliderDate,'Value');
    if isfield(TestData,'Timestamp')
        date_inf = TestData.Timestamp(end);
    else
        date_inf = datenum([2018,1,1]);
    end
    dt = TestData.Timestamp(2) - TestData.Timestamp(1);
    z = get(handles.sliderZoom,'Value');
    z_max = get(handles.sliderZoom,'Max');
    if z == z_max %show all data
        date_end = date_inf;
    elseif z<1
        date_end = min(date_inf,date_1 + 1);
    elseif z<2
        date_end = min(date_inf,date_1 + 7);
    elseif z<3
        date_end = min(date_inf,date_1 + 31);
    elseif z<4
        date_end = min(date_inf,date_1 + 365);
    end
    date = (date_1:dt:date_end)';


    if ~isfield(Plant,'Building')
        Plant.Building = [];
    end
    if ~isfield(Plant,'fluid_loop')
        Plant.fluid_loop = [];
    end
    [forecast,Plant.Generator,Plant.Building] = update_forecast(Plant.Generator,Plant.Building,Plant.fluid_loop,Plant.subNet,Plant.optimoptions,date(2:end),TestData,[]);%% function that creates demand vector with time intervals coresponding to those selected
    n_plot = length(fieldnames(Plant.subNet));
    if isfield(TestData,'Demand')
        outs =  fieldnames(TestData.Demand);
        xi = nnz(TestData.Timestamp<(date_1-1))+1;
        xf = nnz(TestData.Timestamp<(date_end) & TestData.Timestamp>0);
        for i = 1:1:length(outs)
            test_data.Demand.(outs{i}) = TestData.Demand.(outs{i})(xi:xf);
        end
        test_data.Timestamp = TestData.Timestamp(xi:xf);
    else
        test_data = [];
    end
    plot_project([],Plant.Building,Plant.fluid_loop,n_plot,[],forecast,test_data,handles,'Forecast')
end  

%% Network View Tab

%% Market Tab
% --- Executes on button press in BiddingMarket.
function BiddingMarket_Callback(hObject, eventdata, handles)
if get(handles.CommitedMarket,'Value')==1
    a = get(handles.CommitedMarket,'BackgroundColor');
    b = get(handles.BiddingMarket,'BackgroundColor');
    c = get(handles.CommitedMarket,'ForegroundColor');
    d = get(handles.BiddingMarket,'ForegroundColor');
    set(handles.BiddingMarket,'Value',1,'BackgroundColor',a,'ForegroundColor',c);
    set(handles.CommitedMarket,'Value',0,'BackgroundColor',b,'ForegroundColor',d);
else
    set(handles.BiddingMarket,'Value',1); %was already pressed
end


% --- Executes on button press in CommitedMarket.
function CommitedMarket_Callback(hObject, eventdata, handles)
if get(handles.BiddingMarket,'Value')==1
    a = get(handles.BiddingMarket,'BackgroundColor');
    b = get(handles.CommitedMarket,'BackgroundColor');
    c = get(handles.BiddingMarket,'ForegroundColor');
    d = get(handles.CommitedMarket,'ForegroundColor');
    set(handles.CommitedMarket,'Value',1,'BackgroundColor',a,'ForegroundColor',c);
    set(handles.BiddingMarket,'Value',0,'BackgroundColor',b,'ForegroundColor',d);
else
    set(handles.CommitedMarket,'Value',1); %was already pressed
end

function MarketServices(handles)
global Plant
if ~isfield(Plant,'Market') || ~isfield(Plant.Market,'Price')
    Day_Price = [0.2 0.2 0.2 0.2 0.2 0.3 0.3 0.3 0.3 0.3 0.3 0.8 0.8 0.8 0.8 1.3 1.3 0.8 0.8 0.2 0.2 0.2 0.2 0.2]/10;
    Hour_Price = [0.0 0.0 0.0 0.0 0.7 0.8 1.5 1.2 0.4 0.0 0.2 1.0 1.8 2.2 1.7 1.4 1.0 0.7 0.4 0.0 0.0 0.0 0.0 0.0]/10;
    Realtime_Price = [0.03 0.03 0.04 0.04 0.06 0.06 0.05 0.04 0.05 0.1 0.15 0.18 0.16 0.15 0.17 0.12 0.08 0.08 0.07 0.04 0.04 0.03 0.01 0.01];
    Plant.Market.Price = [Day_Price;Hour_Price;Realtime_Price];
    set(handles.MarketTable,'Data',Plant.Market.Price);
end
plotMarginalCapacityCost(handles)

function MarketPopup_Callback(hObject, eventdata, handles)
plotMarginalCapacityCost(handles)

function MarketTable_CellEditCallback(hObject, eventdata, handles)
global Plant
Plant.Market.Price = get(hObject,'Data');

function CapacityTable_CellEditCallback(hObject, eventdata, handles)
global Plant
Plant.Market.Capacity = get(hObject,'Data');

function MarketSlider_Callback(hObject, eventdata, handles)

function MarketSlider_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function Gen_vs_DR_Callback(hObject, eventdata, handles)

function Gen_vs_DR_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function MarketPopup_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% Settings Tab
function Interval_Callback(hObject, eventdata, handles)
global Plant
Plant.optimoptions.Interval = str2double(get(handles.Interval, 'String'));
function Interval_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Resolution_Callback(hObject, eventdata, handles)
global Plant
Plant.optimoptions.Resolution = str2double(get(handles.Resolution, 'String'));
Plant.subNet = [];

function Resolution_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Horizon_Callback(hObject, eventdata, handles)
global Plant
Plant.optimoptions.Horizon = str2double(get(handles.Horizon, 'String'));
Plant.subNet = [];

function Horizon_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function scaletime_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function excessHeat_Callback(hObject, eventdata, handles)
global Plant
Plant.optimoptions.excessHeat = get(hObject, 'Value');
Plant.subNet = [];

% --- Executes on button press in excessCool.
function excessCool_Callback(hObject, eventdata, handles)
global Plant
Plant.optimoptions.excessCool = get(hObject, 'Value');
Plant.subNet = [];

function changingtimesteps_SelectionChangeFcn(hObject, eventdata, handles)
global Plant
switch get(eventdata.NewValue,'Tag')
    case 'constant'
        Plant.optimoptions.tspacing = 'constant';
    case 'manual'
        Plant.optimoptions.tspacing = 'manual';
        prompt = {'Specify a vector of times out of one horizon for each timestep: (ex: .05,.1,.25,.75,1 would result in one timestep at 5hr, 10hr, 25hr, 75hr, and 100hr for a 100 hour horizon)'};
        dlg_title = 'Manual Timesteps';
        num_lines = 1;
        def_ans = {'0.0035,0.0070,0.0105,0.0140,0.0175,0.0210,0.0417,0.0625,0.0833,0.1042,0.125,0.1667,0.2083,0.25,0.3333,0.4167,0.5,0.5833,0.6667,0.75,0.8333,0.9167,1'};
        a = inputdlg(prompt,dlg_title,num_lines,def_ans);
        Plant.optimoptions.manualT = str2double(strsplit(a{:}, ','));
    case 'linear'
        Plant.optimoptions.tspacing = 'constant';
    case 'logarithm'
        Plant.optimoptions.tspacing = 'logarithm';
end
Plant.subNet = [];

% --- Executes when selected object is changed in uipanelSolverMethod.
function ANNradiobutton_Callback(hObject, eventdata, handles)
function NoMixedInteger_Callback(hObject, eventdata, handles)
function MixedInteger_Callback(hObject, eventdata, handles)
function uipanelSolverMethod_SelectionChangedFcn(hObject, eventdata, handles)
global Plant
set(handles.changingtimesteps, 'Visible', 'on');
set(handles.OtherOptions, 'Visible', 'on');
switch get(eventdata.NewValue,'Tag')
    case 'NoMixedInteger'
        Plant.optimoptions.MixedInteger = false;
        Plant.optimoptions.solver = 'quadprog';
    case 'MixedInteger'
        Plant.optimoptions.MixedInteger = true;
        Plant.optimoptions.solver = 'quadprog';
    case 'ANNradiobutton'
        Plant.optimoptions.MixedInteger = false;
        Plant.optimoptions.solver = 'ANN';
    case 'nrel_gurobi'
        set(handles.changingtimesteps, 'Visible', 'off');
        set(handles.OtherOptions, 'Visible', 'off');
        Plant.optimoptions.MixedInteger = false;
        Plant.optimoptions.solver = 'NREL';
        Plant.optiomptions.EnergyPlus = true;
        set(handles.SimBuilding,'value',true)
end

% --- Executes when selected object is changed in SpinningReserve.
function noSpinReserve_Callback(hObject, eventdata, handles)
function SpinReserve_Callback(hObject, eventdata, handles)
function SpinningReserve_SelectionChangeFcn(hObject, eventdata, handles)
global Plant
switch get(eventdata.NewValue,'Tag')
    case 'noSpinReserve'
        Plant.optimoptions.SpinReserve = false;
    case 'SpinReserve'
        Plant.optimoptions.SpinReserve = true;
        Plant.optimoptions.SpinReservePerc = str2double(get(handles.SpinReservePerc, 'String'));
end
Plant.subNet = [];

function SpinReservePerc_Callback(hObject, eventdata, handles)
global Plant
Plant.optimoptions.SpinReservePerc = str2double(get(handles.SpinReservePerc, 'String'));
Plant.subNet = [];

function SpinReservePerc_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in StorageBuff.
function StorageBuff_Callback(hObject, eventdata, handles)
global Plant
val = get(handles.StorageBuff,'Value');
stor = get(handles.StorageBuff,'UserData');
set(handles.editBuffer,'String',num2str(Plant.Generator(stor(val)).VariableStruct.Buffer));
 
% --- Executes during object creation, after setting all properties.
function StorageBuff_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editBuffer_Callback(hObject, eventdata, handles)
global Plant
val = get(handles.StorageBuff,'Value');
stor = get(handles.StorageBuff,'UserData');
Plant.Generator(stor(val)).VariableStruct.Buffer = str2double(get(handles.editBuffer,'String'));
Plant.subNet = [];

function editBuffer_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in SimBuilding.
function SimBuilding_Callback(hObject, eventdata, handles)
global Plant
Plant.optiomptions.EnergyPlus = get(handles.SimBuilding,'value');

function Topt_Callback(hObject, eventdata, handles)
global Plant
Plant.optimoptions.Topt = str2double(get(handles.Topt, 'String'));
Plant.subNet = [];

function Topt_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Tmpc_Callback(hObject, eventdata, handles)
global Plant
Plant.optimoptions.Tmpc = str2double(get(handles.Tmpc, 'String'));
Plant.subNet = [];

function Tmpc_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function scaletime_Callback(hObject, eventdata, handles)
global Plant
%scaletime: the ratio of emulated time to time in the test-data. For example 24 hours of test data can be run in 3 hours with a scaletime of 8. scaletime enlarges any energy storage and slows down the transient response of any generators
Plant.optimoptions.scaletime = str2double(get(handles.scaletime, 'String'));
Plant.subNet = [];

function ForecastMethod_SelectionChangeFcn(hObject, eventdata, handles)
global Plant TestData
if isfield(Plant,'Building') && ~isempty(Plant.Building) && ~isfield(TestData,'Building')
    TestData.Building = load_test_building(Plant.Building,Plant.Network,TestData.Timestamp,TestData.Weather);
end
switch get(eventdata.NewValue,'Tag')
    case 'ARMA'
        Plant.optimoptions.forecast = 'arma';
    case 'ARIMA'
        Plant.optimoptions.forecast = 'ARIMA';
    case 'NeuralNet'
        Plant.optimoptions.forecast = 'NeuralNet';
    case 'Surface'
        Plant.optimoptions.forecast= 'Surface';
    case 'Perfect'
        Plant.optimoptions.forecast = 'Perfect';
    case 'Building'
        Plant.optimoptions.forecast= 'Building';
        TestData = rmfield(TestData,'Building');
end

function editCommandOnOff_Callback(hObject, eventdata, handles)
recordComm(handles)
function editCommandOnOff_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editCommandSet_Callback(hObject, eventdata, handles)
recordComm(handles)
function editCommandSet_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editMeasureSecondary_Callback(hObject, eventdata, handles)
recordComm(handles)
function editMeasureSecondary_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editMeasurePrimary_Callback(hObject, eventdata, handles)
recordComm(handles)
function editMeasurePrimary_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editMeasureInput_Callback(hObject, eventdata, handles)
recordComm(handles)
function editMeasureInput_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editMeasureOnOff_Callback(hObject, eventdata, handles)
recordComm(handles)

function editMeasureOnOff_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function recordComm(handles)
global Plant
nG = length(Plant.Generator);
if isfield(Plant,'Building') && ~isempty(Plant.Building)
    nB = length(Plant.Building);
else
    nB = 0;
end
click = get(handles.CurrentComm,'String');
list = get(handles.uipanelMain1,'UserData');
i = nonzeros((1:(nG+nB))'*strcmp(click,list));
if i<=nG
    Comm.OnOff = str2double(get(handles.editCommandOnOff,'String'));
    Comm.Set = str2double(get(handles.editCommandSet,'String'));

    Measure.OnOff = str2double(get(handles.editMeasureOnOff,'String'));
    Measure.Input = str2double(get(handles.editMeasureInput,'String'));
    Measure.Primary = str2double(get(handles.editMeasurePrimary,'String'));
    Measure.Secondary = str2double(get(handles.editMeasureSecondary,'String'));
    
    Plant.Generator(i).VariableStruct.Comm = Comm;
    Plant.Generator(i).VariableStruct.Measure = Measure;
else
    Measure.Temperature = str2double(get(handles.editMeasurePrimary,'String'));
    Measure.Humidity = str2double(get(handles.editMeasureSecondary,'String'));
    Plant.Building(i-nG).VariableStruct.Measure = Measure;
end

function popupmenuColor_Callback(hObject, eventdata, handles)
handles = guihandles;
list = get(hObject,'String');
color = list{get(hObject,'Value')};
colormap(handles.ResultPlot1,color);
colorVec = colormap(handles.ResultPlot1);
set(handles.figure1,'UserData',colorVec);

list = get(handles.uipanelMain1,'UserData');
colorsPlot = interp1(linspace(0,1,length(colorVec)),colorVec,linspace(0,1,length(list)));
i = 0;
while isfield(handles,strcat('Generator1_',num2str(i+1)))
    i = i+1;
    set(handles.(strcat('Generator1_',num2str(i))),'BackgroundColor',colorsPlot(i,:))
    set(handles.(strcat('Generator3_',num2str(i))),'BackgroundColor',colorsPlot(i,:))
    set(handles.(strcat('Generator5_',num2str(i))),'BackgroundColor',colorsPlot(i,:))
end
nG = i;
while isfield(handles,strcat('Building1_',num2str(i+1-nG)))
    i = i+1;
    set(handles.(strcat('Building1_',num2str(i-nG))),'BackgroundColor',colorsPlot(i,:))
    set(handles.(strcat('Building3_',num2str(i-nG))),'BackgroundColor',colorsPlot(i,:))
    set(handles.(strcat('Building5_',num2str(i-nG))),'BackgroundColor',colorsPlot(i,:))
end

function popupmenuColor_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
