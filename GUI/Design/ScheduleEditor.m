function varargout = ScheduleEditor(varargin)
%SCHEDULEEDITOR MATLAB code file for ScheduleEditor.fig
%      SCHEDULEEDITOR, by itself, creates a new SCHEDULEEDITOR or raises the existing
%      singleton*.
%
%      H = SCHEDULEEDITOR returns the handle to a new SCHEDULEEDITOR or the handle to
%      the existing singleton*.
%
%      SCHEDULEEDITOR('Property','Value',...) creates a new SCHEDULEEDITOR using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to ScheduleEditor_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      SCHEDULEEDITOR('CALLBACK') and SCHEDULEEDITOR('CALLBACK',hObject,...) call the
%      local function named CALLBACK in SCHEDULEEDITOR.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ScheduleEditor

% Last Modified by GUIDE v2.5 07-Mar-2018 23:49:57

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ScheduleEditor_OpeningFcn, ...
                   'gui_OutputFcn',  @ScheduleEditor_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
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


% --- Executes just before ScheduleEditor is made visible.
function ScheduleEditor_OpeningFcn(hObject, eventdata, handles, varargin)
global testSystems SYSINDEX
handles2 = varargin{3};
x = fieldnames(testSystems(SYSINDEX).Building.Schedule);
set(handles.popupSched,'String',x)
units = 'W/ft^2';
if strcmp(varargin{2},'occ')
    val = find(strcmp(x,'occupancy'));
    set(handles.popupSched,'Value',val)
    load = get(handles2.editOccupancy,'String');
    units = 'ppl/ft^2';
elseif strcmp(varargin{2},'light')
    val = find(strcmp(x,'interiorlights'));
    set(handles.popupSched,'Value',val)
    load = get(handles2.editLighting,'String');
else
    val = find(strcmp(x,'equipment'));
    set(handles.popupSched,'Value',val)
    load = get(handles2.editEquipment,'String');
end
set(handles.editLoad,'String',load)
set(handles.loadUnits,'String',units)

set(handles.tableSchedule,'Visible','off')

set(handles.rampfactor,'String','1')
ramp = 1e-4;
data= decideSched(handles); 
plotSched(data,ramp,handles.axesSchedule)

guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = ScheduleEditor_OutputFcn(hObject, eventdata, handles)
building_edit()

function rampfactor_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function popupSched_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editLoad_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function editLoad_Callback(hObject, eventdata, handles)
global testSystems SYSINDEX
% testSystems(SYSINDEX).Building.occupancy = str2double(get(hObject,'String'))/10.76;%convert to m^2
% MainScreen1('calculateBuildingProfile',handles,2)


function ramp_Callback(hObject, eventdata, handles)
global testSystems SYSINDEX
list = get(handles.popupSched,'String');
val = get(handles.popupSched,'Value');
sched = list{val};

if get(hObject,'Value')==1
    set(handles.rampfactor,'Visible','on')
    ramp = str2double(get(handles.rampfactor,'String'));
    testSystems(SYSINDEX).Building.Schedule.(sched).Ramp = ramp;
else
    testSystems(SYSINDEX).Building.Schedule.(sched) = rmfield(testSystems(SYSINDEX).Building.Schedule.(sched),'Ramp');
end
rampfactor_Callback(handles.rampfactor,eventdata,handles)


function step_Callback(hObject, eventdata, handles)
global testSystems SYSINDEX
list = get(handles.popupSched,'String');
val = get(handles.popupSched,'Value');
sched = list{val};
if get(hObject,'Value')==1
    set(handles.rampfactor,'Visible','off')
    testSystems(SYSINDEX).Building.Schedule.(sched) = rmfield(testSystems(SYSINDEX).Building.Schedule.(sched),'Ramp');
end
ramp = 1e-4;
data= decideSched(handles); 
plotSched(data,ramp,handles.axesSchedule)


function rampfactor_Callback(hObject, eventdata, handles)
global testSystems SYSINDEX
list = get(handles.popupSched,'String');
val = get(handles.popupSched,'Value');
sched = list{val};
ramp = str2double(get(handles.rampfactor,'String'));
testSystems(SYSINDEX).Building.Schedule.(sched).Ramp = ramp;
data= decideSched(handles); 
plotSched(data,ramp,handles.axesSchedule)

function onSeason_Callback(hObject, eventdata, handles)

function offSeason_Callback(hObject, eventdata, handles)

function checkEdit_Callback(hObject, eventdata, handles)
if get(handles.checkEdit,'Value')==1
    data = decideSched(handles);
    set(handles.tableSchedule,'Data',data)
    hw = get(handles.tableSchedule,'Extent');
    pos = get(handles.tableSchedule,'Position');
    pos(3:4) = hw(3:4);
    set(handles.tableSchedule,'Visible','on','Position',pos)
else
    set(handles.tableSchedule,'Visible','off')
end

function tableSchedule_CellEditCallback(hObject, eventdata, handles)
global testSystems SYSINDEX
list = get(handles.popupSched,'String');
val = get(handles.popupSched,'Value');
sched = list{val};
if get(handles.weekdaySched,'Value')==1
    type = 'Weekday';
elseif get(handles.saturdaySched,'Value')==1
    type = '.Sat';
else
    type = '.Sun';
end
testSystems(SYSINDEX).Building.Schedule.(sched).(type) = hObject.Data;
if strcmp(get(handles.rampfactor,'Visible'),'on')
    ramp = str2double(get(handles.rampfactor,'String'));
else
    ramp = 1e-4;
end
plotSched(hObject.Data,ramp,handles.axesSchedule)

function popupSched_Callback(hObject, eventdata, handles)
global testSystems SYSINDEX
if strcmp(get(handles.rampfactor,'Visible'),'on')
    ramp = str2double(get(handles.rampfactor,'String'));
else
    ramp = 1e-4;
end
units = 'W/ft^2';
x = get(hObject,'String');
select = get(hObject,'Value');
Value = {'TsetH' 'TsetC' 'daylighting'};
if nnz(strcmp(x(select),Value))
    set(handles.toggleValue,'Value',1)
    toggleValue_Callback(hObject, eventdata, handles)
elseif strcmp(x(select),'exteriorlights_solarcontroled') || strcmp(x(select),'exteriorlights_fixed')
    load=testSystems(SYSINDEX).Building.VariableStruct.ExteriorLights*10.76;
    set(handles.editLoad,'String',num2str(load))
    set(handles.loadUnits,'String',units)
else
    set(handles.toggleFraction,'Value',1)
    toggleFraction_Callback(hObject, eventdata, handles)
	loads = fieldnames(testSystems(SYSINDEX).Building.VariableStruct);
    name = char(x(select));
    name = name(find(~isspace(name)));
    for i=1:length(loads)
        if strcmpi(name,loads(i))
            load = eval(strcat('testSystems(SYSINDEX).Building.VariableStruct.',char(loads(i))))*10.76;
        end
    end 
    set(handles.editLoad,'String',num2str(load))
    if strcmp(x(select),'Occupancy')
        units = 'ppl/ft^2';
    end
    set(handles.loadUnits,'String',units)
end


data= decideSched(handles); 
plotSched(data,ramp,handles.axesSchedule)
checkEdit_Callback(hObject, eventdata, handles)

function toggleFraction_Callback(hObject, eventdata, handles)

set(handles.textLoad,'Visible','on')
set(handles.editLoad,'Visible','on')
set(handles.loadUnits,'Visible','on')



function toggleValue_Callback(hObject, eventdata, handles)

set(handles.textLoad,'Visible','off')
set(handles.editLoad,'Visible','off')
set(handles.loadUnits,'Visible','off')


function weekdaySched_Callback(hObject, eventdata, handles)
if get(hObject,'Value')==1
    set(handles.saturdaySched,'Value',0)
    set(handles.sundaySched,'Value',0)
else
    set(hObject,'Value',1)
end
if strcmp(get(handles.rampfactor,'Visible'),'on')
    ramp = str2double(get(handles.rampfactor,'String'));
else
    ramp = 1e-4;
end
data = decideSched(handles);
plotSched(data,ramp,handles.axesSchedule)
checkEdit_Callback(hObject, eventdata, handles)

function saturdaySched_Callback(hObject, eventdata, handles)
if get(hObject,'Value')==1
    set(handles.weekdaySched,'Value',0)
    set(handles.sundaySched,'Value',0)
else
    set(hObject,'Value',1)
end
if strcmp(get(handles.rampfactor,'Visible'),'on')
    ramp = str2double(get(handles.rampfactor,'String'));
else
    ramp = 1e-4;
end
data= decideSched(handles);
plotSched(data,ramp,handles.axesSchedule)
checkEdit_Callback(hObject, eventdata, handles)

function sundaySched_Callback(hObject, eventdata, handles)
if get(hObject,'Value')==1
    set(handles.saturdaySched,'Value',0)
    set(handles.weekdaySched,'Value',0)
else
    set(hObject,'Value',1)
end
if strcmp(get(handles.rampfactor,'Visible'),'on')
    ramp = str2double(get(handles.rampfactor,'String'));
else
    ramp = 1e-4;
end
data= decideSched(handles); 
plotSched(data,ramp,handles.axesSchedule)
checkEdit_Callback(hObject, eventdata, handles)

function data = decideSched(handles)
global testSystems SYSINDEX
list = get(handles.popupSched,'String');
val = get(handles.popupSched,'Value');
sched = list{val};
if get(handles.weekdaySched,'Value')==1
    type = 'Weekday';
elseif get(handles.saturdaySched,'Value')==1
    type = 'Sat';
else
    type = 'Sun';
end
data = testSystems(SYSINDEX).Building.Schedule.(sched).(type);

function plotSched(data,ramp,axes)
newsched = convertSched(data,ramp);
plot(axes,newsched(:,1),newsched(:,2))


function newsched = convertSched(sched,Ramp)
[m,n] = size(sched);
if m==2
    newsched = sched; %this is constant all day, already made 0 hour and 24 hour
else
    newsched = zeros(2*(m-1),n);
    sched(1,2) = sched(2,2);
    newsched(1,:) = sched(1,:);%hour 0
    newsched(end,:) = sched(m,:);%hour 24
    if Ramp<1e-3
        for i = 2:1:m-1 
            newsched(2*i-2,1) = sched(i);
            newsched(2*i-2,2) = sched(i,2);
            newsched(2*i-1,1) = sched(i)+Ramp;
            newsched(2*i-1,2) = sched(i+1,2);
        end
    else 
        for i = 2:1:m-1 %add points in the middle so it can be properly interpolated
            t_bef = sched(i) - sched(i-1);
            t_aft = sched(i+1)-sched(i);
            newsched(2*i-2,1) = sched(i)-min(Ramp,t_bef/2);
            newsched(2*i-2,2) = sched(i,2);
            newsched(2*i-1,1) = sched(i)+min(Ramp,t_aft/2);
            newsched(2*i-1,2) = sched(i+1,2);
        end
    end
end
