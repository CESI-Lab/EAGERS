function varargout = MainScreen1(varargin)
% MAINSCREEN1 MATLAB code for MainScreen1.fig
%      MAINSCREEN1, by itself, creates a new MAINSCREEN1 or raises the existing
%      singleton*.
%
%      H = MAINSCREEN1 returns the handle to a new MAINSCREEN1 or the handle to
%      the existing singleton*.
%
%      MAINSCREEN1('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MAINSCREEN1.M with the given input arguments.
%
%      MAINSCREEN1('Property','Value',...) creates a new MAINSCREEN1 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MainScreen1_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MainScreen1_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MainScreen1

% Last Modified by GUIDE v2.5 19-Jun-2018 16:55:11

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MainScreen1_OpeningFcn, ...
                   'gui_OutputFcn',  @MainScreen1_OutputFcn, ...
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

% --- Executes just before MainScreen1 is made visible.
function MainScreen1_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MainScreen1 (see VARARGIN)

% Choose default command line output for MainScreen1
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes MainScreen1 wait for user response (see UIRESUME)
% uiwait(handles.figure1);
global Plant SYSINDEX testSystems TestData mainFig
mainFig = gcf;
movegui(gcf,'center');
set(gcf,'Name','Energy Planning Tool 2018.0.2')
%%Plant is necessary because when the saved file opens it is in Plant
%%structure. When returning from EDC it identifies which plant was being
%%simulated
if ~isfield(Plant,'Costs') 
    Costs = [];
else
    Costs = Plant.Costs;
end
Plant.Costs = defaultCosts(Costs,Plant.Generator);
list_proj = {};
if isempty(testSystems)
    SYSINDEX = 1;%index in the list of systems (full Plant structures)
else
    for i = 1:1:length(testSystems)
        list_proj(end+1) = {testSystems(i).Name};
    end
    SYSINDEX = nonzeros((1:length(list_proj))'.*strcmp(Plant.Name,list_proj));
    if isempty(SYSINDEX)
        SYSINDEX = length(testSystems)+1;
    end
end
list_proj(SYSINDEX) = {Plant.Name};
F = fieldnames(Plant);
for j = 1:1:length(F)
    testSystems(SYSINDEX).(F{j}) = Plant.(F{j});
end

%% Main Tabs
%Tags of main tab panels are of form, 'uipanelMain1', 'uipanelMain2', etc.
TabText = {'System Specification';'Cost and Sizing';'Simulation and Analysis';};
for i = 1:length(TabText)
    set(handles.(strcat('MainTab',num2str(i))),'String',TabText{i},'Units','Characters');
    set(handles.(strcat('uipanelMain',num2str(i))),'BorderType','none','Units','Characters')
    if i ==1
        pan1pos = get(handles.uipanelMain1,'Position');
        set(handles.uipanelMain1,'Visible','on')
    else
        set(handles.(strcat('uipanelMain',num2str(i))),'Position',pan1pos,'Visible','off')
    end
end
set(handles.ProjectsPanel,'Units','characters','Position',[110, 41,100,6]);
handles = project_panel_list(list_proj,SYSINDEX,handles);
set(handles.plotting_panel,'Units','characters','Position',[0, 3,215,40]);

handles.axesMain = axes('Units','characters','Position', [13,8,125,24],'NextPlot','add',...
    'Tag', 'axesMain','Parent', handles.plotting_panel,'Visible','on');%primary axes
handles.axesMainR = axes('Units','characters','Position',[13,8,125,24],'NextPlot','add',...
        'Tag', 'axesMainR','Parent', handles.plotting_panel,'color','none',...
        'YAxisLocation','right','Visible','on','xtick',[],'xticklabel',[],'YTick',[],'YTickLabel',[]);
handles.axesCumulative = axes('Units','characters','Position', [155,8,50,24],'NextPlot','add',...
    'color','none','Tag', 'axesCumulative','Parent', handles.plotting_panel,'Visible','on');%histogram axes
colormap(handles.axesMain,'parula');
colormap(handles.axesMainR,'parula');
set(hObject,'UserData',colormap(handles.axesMain));

handles.Optimize = uicontrol('Style', 'pushbutton', 'String', 'Optimize Size',...
    'Units','characters','Position', [105 44 40 3],'BackgroundColor',[0.98,0.5,0.3],...
    'Tag', 'OptimizeSys','FontSize', 12,...
    'Parent', handles.uipanelMain2','Callback','optimize_plant_size','Visible','on');
set(handles.NPC_discount,'String',num2str(testSystems(SYSINDEX).Costs.DiscountRate));

days = max(2,floor(TestData.Timestamp(end) - TestData.Timestamp(1)));
if days<7
    set(handles.sliderZoom1,'Min',0,'Max',1,'Value',0,'SliderStep',[1,1]) %either single day or all data    
elseif days<31
    set(handles.sliderZoom1,'Min',0,'Max',2,'Value',0,'SliderStep',[1/2,1/2]) %either single day, week or all data
elseif days<367
    set(handles.sliderZoom1,'Min',0,'Max',3,'Value',0,'SliderStep',[1/3,1/3]) %either single day, week, month, or all data
else
    set(handles.sliderZoom1,'Min',0,'Max',4,'Value',0,'SliderStep',[1/4,1/4]) %either single day, week, month, year, or all data
end
if days>20
    set(handles.sliderDate1,'Min',0,'Max',days,'Value',0,'SliderStep',[1/days,10/days])
else
    set(handles.sliderDate1,'Min',0,'Max',days,'Value',0,'SliderStep',[1/days,1/days])
end
energy_demands_Callback([], [], handles)
popupmenu_Callback


function popupmenu_Callback
handles = guihandles;
tab = find_active_tab(handles);
if tab == 1
    set(handles.ProjectsPanel,'Visible','on')
    set(handles.plotting_panel,'Visible','off')
    Library_Callback([],[], handles)
    network_representation(handles);
elseif tab == 2
    set(handles.ProjectsPanel,'Visible','off')
    set(handles.plotting_panel,'Visible','off')
    update_costs_table(handles)
elseif tab == 3
    set(handles.ProjectsPanel,'Visible','on')
    set(handles.plotting_panel,'Visible','on','Position',[0,3,215,40]);
    second_plot_vis(handles,'on')
    set(handles.axesMain,'Units','characters','Position',[13,8,125,24]);set(handles.textDate1,'Units','characters','Position',[40,1,8,1.5]);set(handles.sliderDate1,'Units','characters','Position',[48,1,70,1.5]);
end


function load_settings(hObject, eventdata, handles)
global SYSINDEX testSystems
handles = guihandles;
set(handles.uipanelOptimizationOptions,'Visible','on','Parent',get(handles.plotting_panel,'Parent'))
set(handles.ProjectsPanel,'Visible','off')
second_plot_vis(handles,'off')
set(handles.open_settings,'Visible','off')

set(handles.excessHeat, 'value', testSystems(SYSINDEX).optimoptions.excessHeat);
set(handles.excessCool, 'value', testSystems(SYSINDEX).optimoptions.excessCool);
set(handles.NoMixedInteger, 'value', ~testSystems(SYSINDEX).optimoptions.MixedInteger);
set(handles.MixedInteger, 'value', testSystems(SYSINDEX).optimoptions.MixedInteger);

set(handles.SpinReserve, 'value', testSystems(SYSINDEX).optimoptions.SpinReserve);
if testSystems(SYSINDEX).optimoptions.SpinReserve
    set(handles.SpinReservePerc, 'Visible','on','string', testSystems(SYSINDEX).optimoptions.SpinReservePerc);
else set(handles.SpinReservePerc, 'Visible','off');
end

nG = length(testSystems(SYSINDEX).Generator);
str = {};
stor = [];
for i = 1:1:nG
    if ismember(testSystems(SYSINDEX).Generator(i).Type,{'Electric Storage';'Thermal Storage'; 'Hydro Storage';})
        str{end+1} = testSystems(SYSINDEX).Generator(i).Name;
        stor(end+1) = i;
        if ~isfield(testSystems(SYSINDEX).Generator(i).VariableStruct,'Buffer')
            testSystems(SYSINDEX).Generator(i).VariableStruct.Buffer = 20;
        end
    end
end
if ~isempty(str)
    set(handles.StorageBuff, 'Visible', 'on','UserData',stor);
    set(handles.editBuffer, 'Visible', 'on');
    set(handles.textBuffer, 'Visible', 'on');
    set(handles.StorageBuff,'string',str,'value',1);
    set(handles.editBuffer, 'string', testSystems(SYSINDEX).Generator(stor(1)).VariableStruct.Buffer);
else
    set(handles.StorageBuff, 'Visible', 'off');
    set(handles.editBuffer, 'Visible', 'off');
    set(handles.textBuffer, 'Visible', 'off');
end

function close_settings(hObject, eventdata, handles)
handles = guihandles;
set(handles.uipanelOptimizationOptions,'Visible','off')
set(handles.ProjectsPanel,'Visible','on')
second_plot_vis(handles,'on')
set(handles.open_settings,'Visible','on')

% --- Outputs from this function are returned to the command line.
function varargout = MainScreen1_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% Main tabs callback
function mainTab_Callback(hObject, eventdata, handles)
n = get(hObject,'Tag');
n = n(end);
tab = find_active_tab(handles);
tab = num2str(tab);

% CRUCIAL IN NEXT 3 STEPS: m, then n.
% Save color
bColor = get(handles.(strcat('MainTab',tab)),'BackgroundColor');
set(handles.(strcat('MainTab',tab)),'BackgroundColor',max(0,bColor-.1))
bColor = get(handles.(strcat('MainTab',n)),'BackgroundColor');
set(handles.(strcat('MainTab',n)),'BackgroundColor',min(1,bColor+.1))

% Save dimensions
pos = get(handles.(strcat('MainTab',tab)),'Position');
set(handles.(strcat('MainTab',tab)),'Position',[pos(1),pos(2),pos(3),pos(4)-.5])
pos = get(handles.(strcat('MainTab',n)),'Position');
set(handles.(strcat('MainTab',n)),'Position',[pos(1),pos(2),pos(3),pos(4)+.5])

%Change tab font
set(handles.(strcat('MainTab',tab)),'FontWeight','normal','ForegroundColor',[0.501960784313726,0.501960784313726,0.501960784313726])
set(handles.(strcat('MainTab',n)),'FontWeight','bold','ForegroundColor',[0,0,0])

% change visibility
set(handles.(strcat('uipanelMain',tab)),'Visible','off')
set(handles.(strcat('uipanelMain',n)),'Visible','on')
popupmenu_Callback

%% ProjectPanel Functions
function Save_Callback(hObject, eventdata, handles)
global testSystems SYSINDEX Plant Model_dir
Plant = testSystems(SYSINDEX);
[f,p]=uiputfile(fullfile(Model_dir,'Projects','New Project.mat'),'Save Project As...');
if f==0; return; end
Plant.Name=strrep(f,'.mat','');
save([p,f],'Plant')

function Load_Callback(hObject, eventdata, handles)
global Plant Model_dir testSystems SYSINDEX
cd(fullfile(Model_dir,'Projects'))
[fn,pn,~] = uigetfile('*.mat','Load Project File');
load(fullfile(pn,fn));
cd(Model_dir)
if ~isfield(Plant,'Costs') 
    Costs = [];
else
    Costs = Plant.Costs;
end
if isfield(testSystems,'Building') && ~isempty(testSystems(1).Building)
    Plant.Building = testSystems(1).Building;
end
Plant.Costs = defaultCosts(Costs,Plant.Generator);
SYSINDEX = length(testSystems)+1;
F = fieldnames(Plant);
for j = 1:1:length(F)
    testSystems(SYSINDEX).(F{j}) = Plant.(F{j});
end
handles = guihandles;
list_proj={};
for i=1:length(testSystems)
    list_proj(end+1) = {testSystems(i).Name};
end
set(handles.ProjectsPanel,'UserData',list_proj)
handles = project_panel_list(list_proj,SYSINDEX,handles);
System_Callback(handles.(strcat('System_',num2str(SYSINDEX))), eventdata, handles)
popupmenu_Callback
set(handles.popupmenuBaseline,'String',get(handles.ProjectsPanel,'UserData'));

function Copy_Callback(hObject, eventdata, handles)
global Plant testSystems SYSINDEX
Plant = testSystems(SYSINDEX);
Plant.Name = char(inputdlg('New Project Name','Select name',1,{strcat(Plant.Name,'_Alt')}));
testSystems(length(testSystems)+1) = Plant;
handles = guihandles;
list_proj={};
for i=1:length(testSystems)
    list_proj(end+1) = {testSystems(i).Name};
end
set(handles.ProjectsPanel,'UserData',list_proj);
handles = project_panel_list(list_proj,SYSINDEX,handles);
System_Callback(handles.(strcat('System_',num2str(length(testSystems)))), eventdata, handles)
popupmenu_Callback
set(handles.popupmenuBaseline,'String',get(handles.ProjectsPanel,'UserData'));

function System_Callback(hObject, eventdata, handles)
global testSystems SYSINDEX GENINDEX
handles = guihandles;
if isfield(handles,strcat('System_',num2str(SYSINDEX)))
    pos = get(handles.(strcat('System_',num2str(SYSINDEX))),'Position');
    pos(4) = 1.8;
    set(handles.(strcat('System_',num2str(SYSINDEX))),'BackgroundColor',[1 1 .7],'ForegroundColor',[0.5 0.5 0.5],'Position',pos)
end
SYSINDEX = get(hObject,'UserData');
pos = get(handles.(strcat('System_',num2str(SYSINDEX))),'Position');
pos(4) = 2;
set(handles.(strcat('System_',num2str(SYSINDEX))),'BackgroundColor',[1 1 1],'ForegroundColor',[0 0 0],'Position',pos)
pos2 = get(handles.Save,'Position');
pos2(1) = pos(1);
set(handles.Save,'Position',pos2);
GENINDEX = min(length(testSystems(SYSINDEX).Generator), GENINDEX); % prevent index exceeds matrix dimensions error
popupmenu_Callback

function DeleteSys_Callback(hObject, eventdata, handles)
global testSystems SYSINDEX GENINDEX
GENINDEX = 1;
SYSINDEX = get(hObject,'UserData');
val = get(handles.popupmenuBaseline,'Value');
if SYSINDEX<=val
    val = max(1,val-1);
end
if length(testSystems)>SYSINDEX
    testSystems(SYSINDEX:length(testSystems)-1) = testSystems(SYSINDEX+1:length(testSystems));
end
testSystems = testSystems(1:length(testSystems)-1);
handles = guihandles;
list_proj={};
for i=1:length(testSystems)
    list_proj(end+1) = {testSystems(i).Name};
end
SYSINDEX = max(1,SYSINDEX-1);
handles = project_panel_list(list_proj,SYSINDEX,handles);
System_Callback(handles.(strcat('System_',num2str(1))), eventdata, handles)

set(handles.ProjectsPanel,'UserData',list_proj);
set(handles.popupmenuAxes,'Value',1);
set(handles.sliderZoom1,'Visible','on');set(handles.sliderDate1,'Visible','on');set(handles.textDate1,'Visible','on');
set(handles.textDay1,'Visible','on'); set(handles.textAllData1,'Visible','on'); set(handles.textHorizon1,'Visible','on');
set(handles.popupmenuBaseline,'String',get(handles.ProjectsPanel,'UserData'),'Value',val);
update_plots

function PrevSys_Callback(hObject, eventdata, handles)
list = get(handles.ProjectsPanel,'UserData');
gen = length(list);
page = get(handles.PrevSys,'UserData');%current page of the list
if page<2
    set(handles.PrevSys,'Visible','off','UserData',page-1)
else
    set(handles.PrevSys,'Visible','on','UserData',page-1);
end
set(handles.NextSys,'Visible','on','UserData',page-1)
for i = 1:1:12
    if 12*(page-1)+i<=gen
         j = num2str(12*(page-1)+i);
         set(handles.(strcat('System_',j)),'Visible','off');
         set(handles.(strcat('DeleteSys',j)),'Visible','off');
    end
end
for i = 1:1:12
    j = num2str(12*(page-2)+i);
     set(handles.(strcat('System_',j)),'Visible','on');
     set(handles.(strcat('DeleteSys',j)),'Visible','on');
end

function NextSys_Callback(hObject, eventdata, handles)
list = get(handles.ProjectsPanel,'UserData');
gen = length(list);
page = get(handles.PrevSys,'UserData');%current page of the list
if page==ceil(gen/12)-1
    set(handles.NextSys,'Visible','off','UserData',page+1)
else
    set(handles.NextSys,'Visible','on','UserData',page+1);
end
set(handles.PrevSys,'Visible','on','UserData',page+1)
for i = 1:1:12
     j = num2str(12*(page-1)+i);
     set(handles.(strcat('System_',j)),'Visible','off');
     set(handles.(strcat('DeleteSys',j)),'Visible','off');
end
for i = 1:1:12
    j = num2str(12*(page)+i);
    if 12*(page)+i<=gen
         set(handles.(strcat('System_',j)),'Visible','on');
         set(handles.(strcat('DeleteSys',j)),'Visible','on');
    end
end

function pushbuttonEDC_Callback(hObject, eventdata, handles)
global testSystems SYSINDEX Plant mainFig
mainFig = [];  
Plant = testSystems(SYSINDEX);
close
DISPATCH

function tab = find_active_tab(handles)
tab = [];
i = 1;
% Find out which tab is currently selected
while isempty(tab) && isfield(handles,strcat('uipanelMain',num2str(i)))
    if strcmp(get(handles.(strcat('uipanelMain',num2str(i))),'Visible'),'on')
        tab = i;
    else
        i = i+1;
    end
end

%% Tab 1 functions (any changes to the system specification or perfomance features requires the designn days to be re-run)
function saveSystem_Callback(hObject, eventdata, handles)
global GENINDEX SYSINDEX testSystems Model_dir
component = testSystems(SYSINDEX).Generator(GENINDEX);
[f,p]=uiputfile(fullfile(Model_dir,'System Library',component.Name),'Save Plant As...');
if f==0; return; end
save([p,f],'component')

function network_select(hObject, eventdata, handles)
%%show lines connecting the specified network
handles = guihandles;
str = get(handles.network_select,'String');
val = get(handles.network_select,'Value');
net = str{val};
old_val = get(handles.network_select,'UserData');
if iscell(old_val)
    old_net = str{old_val{1}};
else
    old_net = str{old_val};
end
set(handles.(strcat('panel_',old_net)),'Visible','off');
set(handles.(strcat('panel_',net)),'Visible','on');
set(handles.network_select,'UserData',val,'Visible','on');
set(handles.uipanelMain3,'Visible','off')
set(handles.uipanelMain3,'Visible','on')

function node_menu(hObject, eventdata, handles)
global SYSINDEX testSystems GENINDEX
handles = guihandles;
list = get(hObject,'String');
val = get(hObject,'Value');
if val == 1 %zoom into look at node
    str = get(handles.network_select,'String');
    val2 = get(handles.network_select,'Value');
    net = str{val2};
    n_name_fix = strrep(strrep(list{1},' ','_'),',','');
    set(handles.(strcat('node_',n_name_fix)),'Visible','on');
    set(handles.(strcat('panel_',net)),'Visible','off')
else%edit component
    n_g = length(testSystems(SYSINDEX).Generator);
    gen_names = cell(n_g,1);
    for i = 1:1:n_g
        gen_names(i) = {testSystems(SYSINDEX).Generator(i).Name};
    end
    sel_name = list{val};
    GENINDEX = nonzeros(strcmp(sel_name,gen_names));
    if ~isempty(GENINDEX)
        component_details(handles)    
    else
        set(handles.uipanelLibrary,'Visible','off');
        set(handles.building_edit,'Visible','on');
    end
end

function edit_system(hObject, eventdata, handles)
global SYSINDEX testSystems GENINDEX
str = get(hObject,'String');
n_g = length(testSystems(SYSINDEX).Generator);
gen_names = cell(n_g,1);
for i = 1:1:n_g
    gen_names(i) = {testSystems(SYSINDEX).Generator(i).Name};
end
GENINDEX = nonzeros((1:n_g)'.*strcmp(str,gen_names));
if ~isempty(GENINDEX)
    component_details(handles)        
end

function edit_building(hObject, eventdata, handles)
global SYSINDEX testSystems BUILDINDEX
handles = guihandles;
set(handles.uipanelLibrary,'Visible','off');
set(handles.building_edit,'Visible','on');
for i = 1:1:length(testSystems(SYSINDEX).Network)
    n_name_fix = strrep(strrep(testSystems(SYSINDEX).Network(i).name,' ','_'),',','');
    set(handles.(strcat('node_',n_name_fix)),'Visible','off');
end
str = get(handles.network_select,'String');
val = get(handles.network_select,'Value');
net = str{val};
set(handles.(strcat('panel_',net)),'Visible','off')

set(handles.plotting_panel,'Visible','on','Position',[0,0,103,40])
set(handles.axesMain,'Units','characters','Position',[13,8,87,24]);set(handles.textDate1,'Units','characters','Position',[5,1,8,1.5]);set(handles.sliderDate1,'Units','characters','Position',[13,1,70,1.5]);
second_plot_vis(handles,'off')
b_names = {};
for i = 1:1:length(testSystems(SYSINDEX).Building)
    b_names(end+1) = {strrep(strrep(testSystems(SYSINDEX).Building(i).Name,' ','_'),',','')};
end
tag_name = get(hObject,'Tag');
r = strfind(tag_name,'_');
BUILDINDEX = nonzeros((1:length(b_names)).*strcmp({tag_name(r(end)+1:end)},b_names)); 
load_building_parameters(testSystems(SYSINDEX).Building(BUILDINDEX),handles);


function node_view(hObject, eventdata, handles)
global testSystems SYSINDEX
handles = guihandles;
str = get(hObject,'Tag');
str = strrep(str,'_return','');
set(handles.(str),'Visible','off');
set(handles.plotting_panel,'Visible','off','Position',[0,3,215,40])
set(handles.axesMain,'Units','characters','Position',[13,8,125,24]);set(handles.textDate1,'Units','characters','Position',[40,1,8,1.5]);set(handles.sliderDate1,'Units','characters','Position',[48,1,70,1.5]);
if length(testSystems(SYSINDEX).Network) > 1
    network_select([], [], handles)
else
    n_name_fix = strrep(strrep(testSystems(SYSINDEX).Network.name,' ','_'),',','');
    set(handles.(strcat('node_',n_name_fix)),'Visible','on');
end
Library_Callback(hObject, eventdata, handles)

% --- Executes on button press in Library.
function Library_Callback(hObject, eventdata, handles)
set(handles.uipanelLibrary,'Visible','on');
set(handles.uipanelGenSpec,'Visible','off');
set(handles.Library,'Visible','off');
set(handles.saveSystem,'Visible','off');
set(handles.pushbuttonRemove,'Visible','off');

function CompName_Callback(hObject, eventdata, handles)
global testSystems SYSINDEX GENINDEX
testSystems(SYSINDEX).Generator(GENINDEX).Name = get(hObject,'String');
old = char(testSystems(SYSINDEX).Network.Equipment(GENINDEX));
first = strfind(old,'.');
old = old(1:first);
new = cellstr(strcat(old,get(hObject,'String')));
testSystems(SYSINDEX).Network.Equipment(GENINDEX) = new;

function CompName_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function AddSystem(hObject, eventdata,handles)
%Adds a generator or storage system to the project
%% Need to add option to specify node that it is added to
global testSystems SYSINDEX Model_dir GENINDEX
emptyStr = '[None]';
switch get(hObject,'Tag')
    case 'ac_gen'
        folders_add = {'ICE';'mGT'};
    case 'dc_gen'
        folders_add = {'PEM';'SOFC';'MCFC';};
    case 'hydro_dam'
        folders_add = {'HydroelectricDam'};
    case 'electrolyzer'
        folders_add = {'Electrolyzer'};
    case 'renewable'
        folders_add = {'SolarPV';'SolarStirling';'SolarThermal';'Wind';};
    case 'building'
        folders_add = {'Buildings'};
    case 'chiller'
        folders_add = {'Chiller'};
    case 'ab_chiller'
        folders_add = {'AbChiller'};
    case 'air_heater'
        folders_add = {'AirHeater'};
    case 'water_heater'
        folders_add = {'WaterHeater'};
    case 'cool_tower'
        folders_add = {'CoolingTower'};
    case 'cold_storage'
        folders_add = {'ColdStor'};
    case 'hot_storage'
        folders_add = {'HotStor'};
    case 'electric_storage'
        folders_add = {'Battery'};
    case 'hydrogen_storage'
        folders_add = {'HydrogenStorage'};
    case 'Utility'
        folders_add = {'Utility'};
end
list_comp = {};
list_dir = {};
for i = 1:1:length(folders_add)
    files=dir(fullfile(Model_dir,'System Library',folders_add{i},'*.mat'));
    list_comp(end+1:end+length(files))=strrep({files.name},'.mat','');
    list_dir(end+1:end+length(files))=folders_add(i);
end
if isempty(list_comp)
    list_comp{1} = emptyStr;
end
handles = guihandles;
[s,OK] = listdlg('PromptString','Select Model','SelectionMode','single','ListString',list_comp);
if OK && strcmp(get(hObject,'Tag'),'building')
    load(fullfile(Model_dir,'System Library','Buildings',strcat(list_comp{s},'.mat')));
    if ~isfield(testSystems(SYSINDEX),'Building') || isempty(testSystems(SYSINDEX).Building)
        testSystems(SYSINDEX).Building = building;
    else
        if isfield(testSystems(SYSINDEX).Building,'QPform') && ~isfield(building,'QPform')
            building.QPform = [];
            building.Tzone = [];
            building.Twall = [];
            building.Timestamp = [];
        end
        testSystems(SYSINDEX).Building(end+1) = building;
    end
    network_representation(handles);
elseif OK && ~strcmp(list_comp{s},emptyStr)
    componentName = list_comp{s};
    load(fullfile(Model_dir,'System Library',list_dir{s},strcat(componentName,'.mat')));
    if isfield(testSystems(SYSINDEX).Generator(1),'QPform')
        component.QPform = [];
        component.CurrentState = [];
        component.Status = [];
    end
    n_g = length(testSystems(SYSINDEX).Generator);
    gen_names = cell(n_g,1);
    for i = 1:1:n_g
        name = testSystems(SYSINDEX).Generator(i).Name;
        r = strfind(name,'_');
        if ~isempty(r) && all(strcmp(name(r(end)+1:end),{'1';'2';'3';'4';'5';'6';'7';'8';'9';'0';}))
            name = name(1:(r(end)-1));
        end
        gen_names(i) = {name};
    end
    same_name = nonzeros(strcmp(component.Name,gen_names));
    if ~isempty(same_name)
        component.Name = strcat(component.Name,'_',num2str(length(same_name)+1));
    end
    GENINDEX = n_g+1;
    testSystems(SYSINDEX).Generator(GENINDEX) = component;
    type_name = strcat(testSystems(SYSINDEX).Generator(GENINDEX).Type,'.',testSystems(SYSINDEX).Generator(GENINDEX).Name);
    list_nodes = cell(length(testSystems(SYSINDEX).Network)+1,1);
    list_nodes(1) = {'New Node'};
    for i = 1:1:length(testSystems(SYSINDEX).Network)
        list_nodes(i+1) = {testSystems(SYSINDEX).Network(i).name};
    end
    [sel,OK] = listdlg('PromptString','Select Node', 'SelectionMode','single','ListString',list_nodes);
    if OK && sel == 1
        testSystems(SYSINDEX).Network = create_new_node(testSystems(SYSINDEX).Network,type_name);
    elseif OK
        testSystems(SYSINDEX).Network(1).Equipment(end+1) = {type_name};
    end
    network_representation(handles);
    component_details(handles)
end
if ~isfield(testSystems(SYSINDEX),'Costs') 
    Costs = [];
else
    Costs = testSystems(SYSINDEX).Costs;
end
testSystems(SYSINDEX).Costs = defaultCosts(Costs,testSystems(SYSINDEX).Generator);
testSystems(SYSINDEX).Design = [];%empty design day solution

% --- Executes on button press in pushbuttonRemove.
function pushbuttonRemove_Callback(hObject, eventdata, handles)
global testSystems SYSINDEX GENINDEX 
handles = guihandles;
str = strcat(testSystems(SYSINDEX).Generator(GENINDEX).Type,'.',testSystems(SYSINDEX).Generator(GENINDEX).Name);
for n = 1:1:length(testSystems(SYSINDEX).Network)
    testSystems(SYSINDEX).Network(n).Equipment = testSystems(SYSINDEX).Network(n).Equipment(~strcmp(str,testSystems(SYSINDEX).Network(n).Equipment));
end
testSystems(SYSINDEX).Generator = testSystems(SYSINDEX).Generator([1:GENINDEX-1,GENINDEX+1:length(testSystems(SYSINDEX).Generator)]);
testSystems(SYSINDEX).Costs.Equipment = testSystems(SYSINDEX).Costs.Equipment([1:GENINDEX-1,GENINDEX+1:length(testSystems(SYSINDEX).Costs.Equipment)]);
GENINDEX = 1;
network_representation(handles);
Library_Callback([],[], handles);
testSystems(SYSINDEX).Design = [];%empty design day solution

% --- Executes on button press in LoadTestData.
function load_demands(hObject, eventdata, handles)
global TestData testSystems
n_s = length(TestData.Timestamp);
tag_name = get(hObject,'Tag');
%find network abreviation and node # for load
r = strfind(tag_name,'_');
n_name = tag_name(r(1)+1:r(2)-1);
list_node = {};
for i = 1:1:length(testSystems(SYSINDEX).Network)
    list_node(end+1) = {strrep(strrep(testSystems(SYSINDEX).Network(i).name,' ','_'),',','')};
end
node_num = nonzeros((1:length(list_node))'.*strcmp(n_name,list_node));
net = tag_name(r(2)+1:r(3)-1);
nn_list = {'Electrical';'DistrictHeat';'DistrictCool';'Hydro';'DirectCurrent';'CoolingWater';'Transmission1';'Hydrogen';'LiqHydrogen';'Heating2';};
nn_abrev = {'E';'H';'C';'W';'DC';'CW';'E1';'Hy';'LH2';'H2';};
nn_index = nonzeros((1:length(nn_list)).*strcmp(net,nn_list));
abrev = nn_abrev(nn_index);
if ~isempty(testSystems(SYSINDEX).Network(node_num).(net).Load)
    dem_num = testSystems(SYSINDEX).Network(node_num).(net).Load;
else
    if isfield(TestData,'Demand') && isfield(TestData.Demand,abrev)
        dem_num = length(TestData.Demand.(abrev)(1,:))+1;
        testSystems(SYSINDEX).Network(node_num).(net).Load = dem_num;
    else
        dem_num = 1;
    end
end
%load some data and interpolate with timestamp.

TestData.Demand.(net)(:,dem_num) = demand;
for i_ts = 1:1:length(testSystems)
    testSystems(i_ts).Design = [];%empty design day solution
end

function SaveBuilding_Callback(hObject, eventdata, handles)
global testSystems SYSINDEX Model_dir
%savebuilding building type
[f,p]=uiputfile(fullfile(Model_dir,'System Library','Buildings',strcat(testSystems(SYSINDEX).Building.Name,'.mat')),'Save Building As...');
if f==0; return; end
testSystems(SYSINDEX).Building.Name = strrep(f,'.mat','');
building = testSystems(SYSINDEX).Building;
save([p,f],'building')

%% Tab 2 functions
% --- Executes when entered data in editable cell(s) in uitableCosts.
function uitableCosts(hObject, eventdata, handles)
global testSystems SYSINDEX
Data = get(hObject,'Data');
for i = 1:1:length(Data(:,1))
    testSystems(SYSINDEX).Costs.Equipment(i).Name = Data{i,1};
    testSystems(SYSINDEX).Costs.Equipment(i).Cost = Data{i,2};
    testSystems(SYSINDEX).Costs.Equipment(i).OandM = Data{i,3};
    testSystems(SYSINDEX).Costs.Equipment(i).Financed = Data{i,4};
    testSystems(SYSINDEX).Costs.Equipment(i).LoanRate = Data{i,5};
    testSystems(SYSINDEX).Costs.Equipment(i).LoanTerm = Data{i,6};
    if  Data{i,8} ~= testSystems(SYSINDEX).Generator(i).Size
        %update_component_spec
    end
end

function update_costs_table(handles)
global testSystems SYSINDEX
n = length(testSystems(SYSINDEX).Generator);
Costs = cell(n,10);
for i = 1:1:n
    Costs(i,1) = {testSystems(SYSINDEX).Costs.Equipment(i).Name};
    if ~isempty(testSystems(SYSINDEX).Costs.Equipment(i).Cost)
        Costs(i,2) = {testSystems(SYSINDEX).Costs.Equipment(i).Cost};
        Costs(i,3) = {testSystems(SYSINDEX).Costs.Equipment(i).OandM};
        Costs(i,4) = {testSystems(SYSINDEX).Costs.Equipment(i).Financed};
        Costs(i,5) = {testSystems(SYSINDEX).Costs.Equipment(i).LoanRate};
        Costs(i,6) = {testSystems(SYSINDEX).Costs.Equipment(i).LoanTerm};
        Costs(i,7) = {testSystems(SYSINDEX).Generator(i).Size};
        Costs(i,8) = {true};
        Costs(i,9) = {0};
        Costs(i,10) = {2*testSystems(SYSINDEX).Generator(i).Size};
    end
end
set(handles.uitableCosts,'Data',Costs)

%% Tab 3 functions
function popupmenuBaseline_Callback(hObject, eventdata, handles)
%do nothing

function NPC_Years_Callback(hObject, eventdata, handles)
%Do nothing

function NPC_discount_Callback(hObject, eventdata, handles)
global testSystems SYSINDEX
testSystems(SYSINDEX).Costs.DiscountRate = str2double(get(hObject,'String'));

%% Setting functions
function AggressiveOpt_Callback(hObject, eventdata, handles)
if get(handles.AggressiveOpt,'Value')
    set(handles.MedAncilOpt,'Value',0)
    set(handles.RobustAncilOpt,'Value',0)
end

function MedAncilOpt_Callback(hObject, eventdata, handles)
if get(handles.MedAncilOpt,'Value')
    set(handles.AggressiveOpt,'Value',0)
    set(handles.RobustAncilOpt,'Value',0)
end

function RobustAncilOpt_Callback(hObject, eventdata, handles)
if get(handles.RobustAncilOpt,'Value')
    set(handles.AggressiveOpt,'Value',0)
    set(handles.MedAncilOpt,'Value',0)
end

function AutoAncilOpt_Callback(hObject, eventdata, handles)
if get(handles.AutoAncilOpt,'Value')
    set(handles.ManualAncilOpt,'Value',0)
else
    set(handles.ManualAncilOpt,'Value',1)
end

function ManualAncilOpt_Callback(hObject, eventdata, handles)
if get(handles.ManualAncilOpt,'Value')
    set(handles.AutoAncilOpt,'Value',0)
else
    set(handles.AutoAncilOpt,'Value',1)
end

function uipanelOptimizationOptions_SelectionChangeFcn(hObject, eventdata, handles)
global testSystems SYSINDEX
switch get(eventdata.NewValue,'Tag')
    case 'NoMixedInteger'
        testSystems(SYSINDEX).optimoptions.MixedInteger = false;
    case 'MixedInteger'
        testSystems(SYSINDEX).optimoptions.MixedInteger = true;
end

function Horizon_Callback(hObject, eventdata, handles)
global testSystems SYSINDEX
testSystems(SYSINDEX).optimoptions.Horizon = str2double(get(handles.Horizon, 'String'));

function Resolution_Callback(hObject, eventdata, handles)
global testSystems SYSINDEX
testSystems(SYSINDEX).optimoptions.Resolution = str2double(get(handles.Resolution, 'String'));

function CarbonTax_KeyPressFcn(hObject, eventdata, handles)

function DesignDay_Callback(hObject, eventdata, handles)


function StorageBuff_Callback(hObject, eventdata, handles)
global testSystems SYSINDEX
val = get(handles.StorageBuff,'Value');
stor = get(handles.StorageBuff,'UserData');
set(handles.editBuffer,'String',num2str(testSystems(SYSINDEX).Generator(stor(val)).VariableStruct.Buffer));

function editBuffer_Callback(hObject, eventdata, handles)
global testSystems SYSINDEX
val = get(handles.StorageBuff,'Value');
stor = get(handles.StorageBuff,'UserData');
testSystems(SYSINDEX).Generator(stor(val)).VariableStruct.Buffer = str2double(get(handles.editBuffer, 'String'));

function excessHeat_Callback(hObject, eventdata, handles)
global testSystems SYSINDEX
testSystems(SYSINDEX).optimoptions.excessHeat = get(hObject, 'Value');

function excessCool_Callback(hObject, eventdata, handles)
global testSystems SYSINDEX
testSystems(SYSINDEX).optimoptions.excessCool = get(hObject, 'Value');

% --- Executes on button press in SpinReserve.
function SpinReserve_Callback(hObject, eventdata, handles)
global testSystems SYSINDEX
if get(hObject,'Value')
    testSystems(SYSINDEX).optimoptions.SpinReserve = true;
    set(handles.SpinReservePerc,'Visible','on')
else
    testSystems(SYSINDEX).optimoptions.SpinReserve = false;
    set(handles.SpinReservePerc,'Visible','off')
end

function SpinReservePerc_Callback(hObject, eventdata, handles)
global testSystems SYSINDEX
testSystems(SYSINDEX).optimoptions.SpinReservePerc = str2double(get(handles.SpinReservePerc, 'String'));


% --- Executes on button press in cost_analysis.
function cost_analysis_Callback(hObject, eventdata, handles)
%update list of plotting options, then 
list = {'Monthly Costs';};
handles = guihandles;
second_plot_vis(handles,'on')
set(handles.popupmenuAxes,'String',list,'Value',1);
set(handles.NPC_discount,'Visible','on');set(handles.text_discount_rate,'Visible','on');set(handles.text_perc,'Visible','on');
set(handles.NPC_Years,'Visible','on');set(handles.text_cost_horizon,'Visible','on');set(handles.text_years,'Visible','on');
if ischar(get(handles.axesMainR,'YLabel'))
    set(handles.axesMainR,'xtick',[],'xticklabel',[],'YTick',[],'YTickLabel',[],'Ylabel',[]);
else
    set(handles.axesMainR,'xtick',[],'xticklabel',[],'YTick',[],'YTickLabel',[]);
end
update_plots

% --- Executes on button press in equipment_dispatch.
function equipment_dispatch_Callback(hObject, eventdata, handles)
global testSystems SYSINDEX
handles = guihandles;
set(handles.sliderZoom1,'Value',0);
second_plot_vis(handles,'off')
% list = {};
% for i = 1:1:length(testSystems(SYSINDEX).Network)
%     list(end+1) = {testSystems(SYSINDEX).Network(i).name};
% end
network_names = fieldnames(testSystems(SYSINDEX).Network);
network_names = network_names(~strcmp('name',network_names));
network_names = network_names(~strcmp('Equipment',network_names));
network_names = network_names(~strcmp('Location',network_names));
list = network_names;
if isfield(testSystems(SYSINDEX),'Building')
    for i = 1:1:length(testSystems(SYSINDEX).Building)
        list(end+1) = {testSystems(SYSINDEX).Building(i).Name};
    end
end
set(handles.popupmenuAxes,'String',list,'Value',1);
set(handles.NPC_discount,'Visible','off');set(handles.text_discount_rate,'Visible','off');set(handles.text_perc,'Visible','off');
set(handles.NPC_Years,'Visible','off');set(handles.text_cost_horizon,'Visible','off');set(handles.text_years,'Visible','off');
update_plots

% --- Executes on button press in energy_demands.
function energy_demands_Callback(hObject, eventdata, handles)
global TestData testSystems SYSINDEX
handles = guihandles;
set(handles.open_settings,'Visible','on')
second_plot_vis(handles,'on')
list = {};
if isfield(TestData,'Demand')
    if isfield(TestData.Demand,'E')
        list(end+1) = {'Electrical Demand'};
    end
    if isfield(TestData.Demand,'H')
        list(end+1) = {'Heating Demand'};
    end
    if isfield(TestData.Demand,'C')
        list(end+1) = {'Cooling Demand'};
    end
    if isfield(TestData.Demand,'W')
        list(end+1) = {'Water Demand'};
    end
end
if isfield(testSystems(SYSINDEX),'Building')
    list(end+1:end+2) = {'InternalGains';'NonHVACelectric';};
    TestData.Building = load_test_building(testSystems(SYSINDEX).Building,testSystems(SYSINDEX).Network,TestData.Timestamp,TestData.Weather);
end
set(handles.popupmenuAxes,'String',list,'Value',1);
set(handles.NPC_discount,'Visible','off');set(handles.text_discount_rate,'Visible','off');set(handles.text_perc,'Visible','off');
set(handles.NPC_Years,'Visible','off');set(handles.text_cost_horizon,'Visible','off');set(handles.text_years,'Visible','off');
if ischar(get(handles.axesMainR,'YLabel'))
    set(handles.axesMainR,'xtick',[],'xticklabel',[],'YTick',[],'YTickLabel',[],'Ylabel',[]);
else
    set(handles.axesMainR,'xtick',[],'xticklabel',[],'YTick',[],'YTickLabel',[]);
end
update_plots


function second_plot_vis(handles,vis)
a = get(handles.axesCumulative,'Children');
set(handles.axesCumulative,'Visible',vis)
for i = 1:1:length(a)
    if iscell(a)
        set(a{i},'Visible',vis);
    else
        set(a(i),'Visible',vis)
    end
end
