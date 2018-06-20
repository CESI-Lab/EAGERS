function building_edit(varargin)
global testSystems SYSINDEX BUILDINDEX TestData
if ~isempty(varargin)
    build = testSystems(SYSINDEX).Building(BUILDINDEX);
    if ischar(varargin{1})
        f= str2func(varargin{1});
        testSystems(SYSINDEX).Building(BUILDINDEX) = f(build,varargin{2:end});
    else
        f = str2func(get(varargin{1},'Call'));
        testSystems(SYSINDEX).Building(BUILDINDEX) = f(build,varargin);
    end
end
TestData.Building = load_test_building(testSystems(SYSINDEX).Building,testSystems(SYSINDEX).Network,TestData.Timestamp,TestData.Weather);
update_plots
testSystems(SYSINDEX).Design = [];%empty design day solution
% end%ends function building_edit

%% Buiding specification
function building = building_list(build,hObject, eventdata, handles)
%load the selected building
global Model_dir 
list = get(handles.building_list,'string');
sel = get(handles.building_list,'Value');
load(fullfile(Model_dir,'System Library','Buildings',list{sel}));
load_building_parameters(handles)

function build = climate_list(build,hObject, eventdata, handles)
%load selected weather profile
global Model_dir TestData
list = get(handles.climate_list,'string');
sel = get(handles.climate_list,'Value');
load(fullfile(Model_dir,'System Library','Weather',list{sel}));
TestData.Weather = interpolate_weather(weather,TestData.Timestamp);

function build = editName(build,hObject, eventdata, handles)
build.Name = get(hObject,'String');

function build = editArea(build,hObject, eventdata, handles)
build.Area = str2double(get(hObject,'String'))/10.76;%convert to m^2
% volume of treated air space (assume height of 3m)
build.Volume = build.Area*3; % m^3

function build = editOccupancy(build,hObject, eventdata, handles)
build.VariableStruct.occupancy = str2double(get(hObject,'String'))/10.76;%convert to m^2

function build = editLighting(build,hObject, eventdata, handles)
build.VariableStruct.lighting = str2double(get(hObject,'String'))/10.76;%convert to m^2

function build = editEquipment(build,hObject, eventdata, handles)
build.VariableStruct.equipment = str2double(get(hObject,'String'))/10.76;%convert to m^2

function build = editComfort(build,hObject, eventdata, handles)
build.VariableStruct.Comfort = str2double(get(hObject,'String'))*5/9;%convert to Celcius

function build = editRvalue(build,hObject, eventdata, handles)
%aproximate R-values: Windows = 3, doors = 5, walls = .03, floor/ceiling 50
% U = 1/R-value
% Q= U*A*deltaT, resistance = 1/U
% a value of U = 0.2 BTU/hr-F-ft^2 converts to 1.1352e-3 kJ/s-K-m^2 (0.2*1055/3600*9/5*10.76/1000)  10.76ft^2 per m^2, 1055J/BTU, 9/5F per K
% inverting gives Resistance = 880.92 m^2*K/kW
build.VariableStruct.WallRvalue = str2double(get(hObject,'String'))*3600/1055*5/9*1000/10.76;% convert hr-F-ft^2/BTU to m^2*K/kW
build.VariableStruct.RoofRvalue = str2double(get(hObject,'String'))*3600/1055*5/9*1000/10.76;% convert hr-F-ft^2/BTU to m^2*K/kW

function build = editAirChange(build,hObject, eventdata, handles)
build.VariableStruct.AirChangePerHr = str2double(get(hObject,'String'));

function build = editDewPoint(build,hObject, eventdata, handles)
build.VariableStruct.DPset = (str2double(get(hObject,'String'))-32)*5/9;

function build = editColdSupply(build,hObject, eventdata, handles)
build.VariableStruct.ColdAirSet = (str2double(get(hObject,'String'))-32)*5/9;

function build = editFanPower(build,hObject, eventdata, handles)
build.VariableStruct.Fan_Power = str2double(get(hObject,'String'));

function build = editWindowWall(build,hObject, eventdata, handles)
build.VariableStruct.WindowWallRatio = str2double(get(hObject,'String'));

function build = editUvalue(build,hObject, eventdata, handles)
build.VariableStruct.WindowUvalue = str2double(get(hObject,'String'));

function build = editDailyVariability(build,hObject, eventdata, handles)
build.VariableStruct.DailyVariability = str2double(get(hObject,'String'));

function build = editHourlyVariability(build,hObject, eventdata, handles)
build.VariableStruct.HourlyVariability = str2double(get(hObject,'String'));

function build = editSubHourlyVariability(hObject, eventdata, handles)
build.VariableStruct.SubHourlyVariability = str2double(get(hObject,'String'));