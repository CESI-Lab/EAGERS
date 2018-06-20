function varargout = WelcomeScreen(varargin)
% WELCOMESCREEN MATLAB code for WelcomeScreen.fig
%      WELCOMESCREEN, by itself, creates a new WELCOMESCREEN or raises the existing
%      singleton*.
%
%      H = WELCOMESCREEN returns the handle to a new WELCOMESCREEN or the handle to
%      the existing singleton*.
%
%      WELCOMESCREEN('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in WELCOMESCREEN.M with the given input arguments.
%
%      WELCOMESCREEN('Property','Value',...) creates a new WELCOMESCREEN or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before WelcomeScreen_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to WelcomeScreen_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help WelcomeScreen

% Last Modified by GUIDE v2.5 18-Feb-2017 17:56:45

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @WelcomeScreen_OpeningFcn, ...
                   'gui_OutputFcn',  @WelcomeScreen_OutputFcn, ...
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


% --- Executes just before WelcomeScreen is made visible.
function WelcomeScreen_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to WelcomeScreen (see VARARGIN)

% Choose default command line output for WelcomeScreen
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes WelcomeScreen wait for user response (see UIRESUME)
% uiwait(handles.figure1);
global testSystems mainFig TestData 
if ishandle(mainFig)
    close(mainFig)
end
mainFig = gcf;
testSystems =[];
TestData = [];
handles.output = hObject;
guidata(hObject, handles);
movegui(gcf,'center');
model_dir=strrep(which('WelcomeScreen.m'),fullfile('GUI','WelcomeScreen.m'),'');
files=dir(fullfile(model_dir, 'Projects','*.mat'));
list=strrep({files.name},'.mat','');
set(handles.ProjectList,'string',list,'value',1)

files = dir(fullfile(model_dir, 'Model Library','*.mat'));
list=strrep({files.name},'.mat','');
set(handles.popupmenuSTRIDES,'string',list,'value',1)

set(gcf,'Name','EAGERS 2018.2.3')


% --- Outputs from this function are returned to the command line.
function varargout = WelcomeScreen_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in ProjectList.
function ProjectList_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function ProjectList_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbuttonDesign.
function pushbuttonDesign_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonDesign (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Plant mainFig
% Load file that was selected from the popupmenu
model_dir=strrep(which('WelcomeScreen.m'),fullfile('GUI','WelcomeScreen.m'),'');
projList = get(handles.ProjectList,'String');
projName = projList{get(handles.ProjectList,'Value')};
load(fullfile(model_dir,'Projects',projName));
load_test_data
mainFig = [];
close
MainScreen1

% --- Executes on button press in Open.
function Open_Callback(hObject, eventdata, handles)
global Plant mainFig
model_dir=strrep(which('WelcomeScreen.m'),fullfile('GUI','WelcomeScreen.m'),'');
list=get(handles.ProjectList,'string');
plantSel = list{get(handles.ProjectList,'value')};
load(fullfile(model_dir,'Projects',plantSel))
load_test_data
mainFig = [];
close
%open new GUI
DISPATCH


% --- Executes on button press in pushbuttonSTRIDES.
function pushbuttonSTRIDES_Callback(hObject, eventdata, handles)
global mainFig
mainFig = [];
close
STRIDES


% --- Executes on selection change in popupmenuSTRIDES.
function popupmenuSTRIDES_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function popupmenuSTRIDES_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
