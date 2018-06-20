%%move to seperate function to handle building stuff
function load_building_parameters(build,handles)
global Model_dir TestData
files = dir(fullfile(Model_dir, 'System Library','Buildings','*.mat'));
listB=strrep({files.name},'.mat','');
set(handles.building_list,'string',listB,'value',1);

files = dir(fullfile(Model_dir, 'System Library','Weather','*.mat'));
listW=strrep({files.name},'.mat','');
set(handles.climate_list,'string',listW,'value',1);
I = find(strcmp(build.Name,listB));
if isempty(I)
    listB = {build.Name;listB};
    I = 1;
end
set(handles.building_list,'string',listB);
set(handles.building_list,'value',I);

I = find(strcmp(TestData.Weather.Name,listW));
if isempty(I)
    listW(end+1) = {TestData.Weather.Name};
    I = length(listW);
end
set(handles.climate_list,'string',listW);
set(handles.climate_list,'value',I);


set(handles.build_name,'String',build.Name)
set(handles.editArea,'String',num2str(build.Area*10.76))
set(handles.editWindowWall,'String',num2str(build.VariableStruct.WindowWallRatio))
set(handles.editOccupancy,'String',num2str(build.VariableStruct.occupancy*10.76))
set(handles.editLighting,'String',num2str(build.VariableStruct.InteriorLights*10.76))
set(handles.editEquipment,'String',num2str(build.VariableStruct.equipment*10.76))
set(handles.editComfort,'String',num2str(build.VariableStruct.Comfort*9/5))
set(handles.editRvalue,'String',num2str(build.VariableStruct.WallRvalue*1055/3600*9/5*10.76/1000))
set(handles.editUvalue,'String',num2str(build.VariableStruct.WindowUvalue))
set(handles.editAirChange,'String',num2str(build.VariableStruct.AirChangePerHr))
set(handles.editDewPoint,'String',num2str(build.VariableStruct.DPset*9/5+32))
set(handles.editColdSupply,'String',num2str(build.VariableStruct.ColdAirSet*9/5+32))
set(handles.editFanPower,'String',num2str(build.VariableStruct.FanPower))

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
update_plots