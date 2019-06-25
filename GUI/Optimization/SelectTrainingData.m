function varargout = SelectTrainingData(varargin)
% SELECTTRAININGDATA MATLAB code for SelectTrainingData.fig
%      SELECTTRAININGDATA, by itself, creates a new SELECTTRAININGDATA or raises the existing
%      singleton*.
%
%      H = SELECTTRAININGDATA returns the handle to a new SELECTTRAININGDATA or the handle to
%      the existing singleton*.
%
%      SELECTTRAININGDATA('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SELECTTRAININGDATA.M with the given input arguments.
%
%      SELECTTRAININGDATA('Property','Value',...) creates a new SELECTTRAININGDATA or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SelectTrainingData_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SelectTrainingData_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SelectTrainingData

% Last Modified by GUIDE v2.5 22-Feb-2018 17:19:13

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SelectTrainingData_OpeningFcn, ...
                   'gui_OutputFcn',  @SelectTrainingData_OutputFcn, ...
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


% --- Executes just before SelectTrainingData is made visible.
function SelectTrainingData_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SelectTrainingData (see VARARGIN)

% Choose default command line output for SelectTrainingData
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

global Model_dir
files = dir(fullfile(Model_dir,'GUI\Optimization\Results','*.mat'));
list = strrep({files.name},'.mat','');
set(handles.TrainingData,'string',list,'value',1)

ANNs = dir(fullfile(Model_dir,'Optimization\NeuralNetwork','*trained.mat'));
ANNlist = strrep({ANNs.name},'.mat','');
set(handles.pretrainedANNList,'string',ANNlist,'value',1)

%UIWAIT makes SelectTrainingData wait for user response (see UIRESUME)
uiwait(handles.figure1);
%TrainFile = get(handles.TrainingData.UserData);
%set(handles.output.UserData, TrainFile);



% --- Outputs from this function are returned to the command line.
function varargout = SelectTrainingData_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
filename = handles.TrainingData.UserData;
varargout{1} = filename;
close


% --- Executes on selection change in TrainingData.
function TrainingData_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function TrainingData_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Okay_Train.
function Okay_Train_Callback(hObject, eventdata, handles)
% hObject    handle to Okay_Train (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
list = get(handles.TrainingData,'string');
file = list{get(handles.TrainingData,'value')};
handles.TrainingData.UserData = strcat('GUI\Optimization\Results\',file,'.mat');
uiresume
%load file that was selected from the popupmenu


% --- Executes on button press in RunNewTraining.
function RunNewTraining_Callback(hObject, eventdata, handles)
% hObject    handle to RunNewTraining (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.TrainingData.UserData = [];
uiresume


% --- Executes on selection change in pretrainedANNList.
function pretrainedANNList_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function pretrainedANNList_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ANNbutton.
function ANNbutton_Callback(hObject, eventdata, handles)
% hObject    handle to ANNbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
list = get(handles.pretrainedANNList,'string');
file = list{get(handles.pretrainedANNList,'value')};
handles.TrainingData.UserData = strcat('Optimization\NeuralNetwork\',file,'.mat');
uiresume
%load ANN that was selected from popupmenu
