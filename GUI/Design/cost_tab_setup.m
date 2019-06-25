function handles = cost_tab_setup(handles)
%% Not USED

%make items on cost tab
list = {'Name';'Size';'Cost';'CostkW';'O_M';'Financed';'Loan';'LoanTerm';};
listText = {'Name';'Size';'Total Cost';'Normalized Cost';'Operations & Maintenance';'Percent Financed';'Interest Rate';'Loan period';};
listUnits = {'';'kW';'$';'$/kW';'$/month';'%';'%';'years';};

quote = '''';
for i=1:1:length(list)
    handles.(strcat('Cost_',list{i},'_text')) = uicontrol('Style', 'text', 'String', listText{i},...
        'Units','characters','Position', [50 49-5*i 35 1.5],'HorizontalAlignment','left',...
        'Tag', strcat('Cost_',list{i},'_text'),'FontSize', 12,'FontWeight','bold',...
        'Parent', handles.uipanelMain2','Visible','off');
    callback = strcat('@(hObject,eventdata)MainScreen1(',quote,'EditGenCost_Callback',quote,',hObject,eventdata,guidata(hObject))');
    handles.(strcat('Cost_',list{i},'_edit')) = uicontrol('Style', 'edit', 'String', '',...
        'Units','characters','Position', [50 47-5*i 20 2],'BackgroundColor',[1,1,1],...
        'Tag', strcat('Cost_',list{i},'_edit'),'FontSize', 10,...
        'Parent', handles.uipanelMain2','UserData',i,'Enable','on',...
        'Callback',eval(callback),'Visible','off');
    handles.(strcat('Cost_',list{i},'_unit')) = uicontrol('Style', 'text', 'String', listUnits{i},...
        'Units','characters','Position', [71 47.25-5*i 25 1.5],'HorizontalAlignment','left',...
        'Tag', strcat('Cost_',list{i},'_unit'),'FontSize', 10,...
        'Parent', handles.uipanelMain2','Visible','off');
end
handles.min_size = uicontrol('Style', 'edit', 'String', '0',...
    'Units','characters','Position', [80 33 15 2],'BackgroundColor',[1,1,1],...
    'Tag', 'min_size','FontSize', 10,'Parent', handles.uipanelMain2','Visible','off');
handles.min_size_text = uicontrol('Style', 'text', 'String', 'Min Size',...
    'Units','characters','Position', [80 31 15 1.5],'HorizontalAlignment','left',...
    'Tag','min_size_text','FontSize', 10,'Parent', handles.uipanelMain2','Visible','off');
handles.max_size = uicontrol('Style', 'edit', 'String', '0',...
    'Units','characters','Position', [100 33 15 2],'BackgroundColor',[1,1,1],...
    'Tag', 'max_size','FontSize', 10,'Parent', handles.uipanelMain2','Visible','off');
handles.max_size_text = uicontrol('Style', 'text', 'String', 'Max Size',...
    'Units','characters','Position', [100 31 15 1.5],'HorizontalAlignment','left',...
    'Tag','max_size_text','FontSize', 10,'Parent', handles.uipanelMain2','Visible','off');

callback = strcat('@(hObject,eventdata)MainScreen1(',quote,'return_to_system',quote,',hObject,eventdata,guidata(hObject))');
handles.ReturnToSys = uicontrol('Style', 'pushbutton', 'String', 'Return to system',...
    'Units','characters','Position', [50 2 35 3],'BackgroundColor',[0.98,0.5,0.3],...
    'Tag', 'ReturnToSys','FontSize', 12,...
    'Parent', handles.uipanelMain2','Callback',eval(callback),'Visible','off');

handles.OptimizeSys = uicontrol('Style', 'pushbutton', 'String', 'Optimize System Size',...
    'Units','characters','Position', [80 21 35 3],'BackgroundColor',[0.98,0.5,0.3],...
    'Tag', 'OptimizeSys','FontSize', 11,...
    'Parent', handles.uipanelMain2','Callback','optimize_plant_size','Visible','on');
handles.OptimizeGen = uicontrol('Style', 'pushbutton', 'String', 'Optimize Component Size',...
    'Units','characters','Position', [80 36 35 3],'BackgroundColor',[1,0,0],...
    'Tag', 'OptimizeGen','FontSize', 11,...
    'Parent', handles.uipanelMain2','Callback','optimize_component_size','Visible','off');



function PrevGen_Callback(hObject, eventdata, handles)
handles = guihandles;
list = get(handles.uipanelMain2,'UserData');
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
         set(handles.(strcat('Equipment_',j)),'Visible','off');
    end
end
for i = 1:1:12
    j = num2str(12*(page-2)+i);
     set(handles.(strcat('Equipment_',j)),'Visible','on');
end

function NextGen_Callback(hObject, eventdata, handles)
handles = guihandles;
list = get(handles.uipanelMain2,'UserData');
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
     set(handles.(strcat('Equipment_',j)),'Visible','off');
end
for i = 1:1:12
    j = num2str(12*(page)+i);
    if 12*(page)+i<=gen
         set(handles.(strcat('Equipment_',j)),'Visible','on');
    end
end

function SetGenCost_Callback(hObject, eventdata, handles)
global testSystems SYSINDEX GENINDEX
if ~isempty(hObject)
    GENINDEX = get(hObject,'UserData');
end
handles = guihandles;
set(handles.OptimizeSys,'Visible','off');
set(handles.OptimizeGen,'Visible','on');
set(handles.min_size,'Visible','on');
set(handles.min_size_text,'Visible','on');
set(handles.max_size,'Visible','on','String',num2str(2*testSystems(SYSINDEX).Generator(GENINDEX).Size));
set(handles.max_size_text,'Visible','on');
set(handles.ReturnToSys,'Visible','on');
set(handles.uitableCosts,'Visible','off');
set(handles.range_table,'Visible','off');

list = {'Name';'Size';'Cost';'CostkW';'O_M';'Financed';'Loan';'LoanTerm';};
listValues = cell(8,1);
listValues(1) = {testSystems(SYSINDEX).Costs.Equipment(GENINDEX).Name};
listValues(2) = {num2str(testSystems(SYSINDEX).Generator(GENINDEX).Size)};
listValues(3) = {num2str(testSystems(SYSINDEX).Costs.Equipment(GENINDEX).Cost)};
listValues(4) = {num2str(testSystems(SYSINDEX).Costs.Equipment(GENINDEX).Cost/testSystems(SYSINDEX).Generator(GENINDEX).Size)};
listValues(5) = {num2str(testSystems(SYSINDEX).Costs.Equipment(GENINDEX).OandM)};
listValues(6) = {num2str(testSystems(SYSINDEX).Costs.Equipment(GENINDEX).Financed)};
listValues(7) = {num2str(testSystems(SYSINDEX).Costs.Equipment(GENINDEX).LoanRate)};
listValues(8) = {num2str(testSystems(SYSINDEX).Costs.Equipment(GENINDEX).LoanTerm)};
for i = 1:1:length(list(:,1))
    set(handles.(strcat('Cost_',list{i},'_text')),'Visible','on');
    set(handles.(strcat('Cost_',list{i},'_edit')), 'String', listValues{i},'Visible','on');
    set(handles.(strcat('Cost_',list{i},'_unit')),'Visible','on');
end


function return_to_system(hObject, eventdata, handles)
handles = guihandles;
set(handles.OptimizeSys,'Visible','on');
set(handles.OptimizeGen,'Visible','off');
set(handles.min_size,'Visible','off');
set(handles.min_size_text,'Visible','off');
set(handles.max_size,'Visible','off');
set(handles.max_size_text,'Visible','off');
set(handles.ReturnToSys,'Visible','off');
set(handles.uitableCosts,'Visible','on');
set(handles.range_table,'Visible','on');
list = {'Name';'Size';'Cost';'CostkW';'O_M';'Financed';'Loan';'LoanTerm';};
for i = 1:1:length(list(:,1))
    set(handles.(strcat('Cost_',list{i},'_text')),'Visible','off');
    set(handles.(strcat('Cost_',list{i},'_edit')),'Visible','off');
    set(handles.(strcat('Cost_',list{i},'_unit')),'Visible','off');
end


function EditGenCost_Callback(hObject, eventdata, handles)
global testSystems SYSINDEX GENINDEX
handles = guihandles;
i = get(hObject,'UserData');
if i == 1
    testSystems(SYSINDEX).Costs.Equipment(GENINDEX).Name = get(hObject,'String');
elseif i == 2
    testSystems(SYSINDEX).Generator(GENINDEX) = update_component_spec(testSystems(SYSINDEX).Generator(GENINDEX),'UB',str2double(get(hObject,'String')));
elseif i == 3
    testSystems(SYSINDEX).Costs.Equipment(GENINDEX).Cost = str2double(get(hObject,'String'));
    set(handles.Cost_CostkW_edit,'String',num2str(testSystems(SYSINDEX).Costs.Equipment(GENINDEX).Cost/testSystems(SYSINDEX).Generator(GENINDEX).Size));
elseif i == 4
    testSystems(SYSINDEX).Costs.Equipment(GENINDEX).Cost = str2double(get(hObject,'String'))*testSystems(SYSINDEX).Generator(GENINDEX).Size;
    set(handles.Cost_Cost_edit,'String',num2str(testSystems(SYSINDEX).Costs.Equipment(GENINDEX).Cost)); 
elseif i == 5
    testSystems(SYSINDEX).Costs.Equipment(GENINDEX).OandM = str2double(get(hObject,'String'));
elseif i == 6
    testSystems(SYSINDEX).Costs.Equipment(GENINDEX).Financed = str2double(get(hObject,'String'));
elseif i == 7
    testSystems(SYSINDEX).Costs.Equipment(GENINDEX).LoanRate = str2double(get(hObject,'String'));
elseif i == 8
    testSystems(SYSINDEX).Costs.Equipment(GENINDEX).LoanTerm = str2double(get(hObject,'String'));
end
updateSummaryTable(handles)