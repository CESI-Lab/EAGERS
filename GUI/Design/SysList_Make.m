function handles = SysList_Make(handles)
%% Create the buttons, check boxes and delete icons for the projects being compared
global testSystems SYSINDEX
quote='''';
list={};
for i=1:length(testSystems)
    list(end+1) = {testSystems(i).Name};
end
%Delete any old buttons
% handlesAll = guihandles;
% handlesAll.axesMain = handles.axesMain;
% handlesAll.axesCumulative = handles.axesCumulative;
% handles = handlesAll; %workaround to pass both axes handles and showSys handles
num = num2str(length(list)+1);
if isfield(handles,strcat('System_',num))
    delete(handles.(strcat('System_',num)));
%     delete(handles.(strcat('showSys',num)));
%     delete(handles.(strcat('textSys',num)));
    delete(handles.(strcat('DeleteSys',num)));
    handles = rmfield(handles,strcat('System_',num));
%     handles = rmfield(handles,strcat('DeleteSys_',num));
end

%make new buttons
nSys = length(list);
for i=1:1:nSys 
    num = num2str(i);
    if ~isfield(handles,strcat('System_',num))
        curtab = floor((i-1)/3)+1;
        prev = 3*(curtab-1);
        xpos = -20 + 24*(i-prev);
        if curtab==1
            vis = 'on';
        else
            vis = 'off';
        end
        if i == SYSINDEX
            foregroundC = [0 0 0];
            backgroundC = [1 1 1];
            height = 2;
        else
            foregroundC = [0.5 0.5 0.5];
            backgroundC = [1 1 .7];
            height = 1.8;
        end
        handles.(strcat('System_',num)) = uicontrol('Style', 'pushbutton', 'String', list{i},...
        'Units','characters',...
        'Position', [xpos 2.5 20 height],...
        'Tag', strcat('System_',num),...
        'FontSize', 10,...
        'Parent', handles.ProjectsPanel,...
        'Callback',eval(strcat('@(hObject,eventdata)MainScreen1(',quote,'System_Callback',quote,',hObject,eventdata,guidata(hObject))')),...
        'Visible',vis,...
        'UserData',i,...
        'ForegroundColor', foregroundC,'BackgroundColor',backgroundC);
%         %% checkbox to show
%         handles.(strcat('showSys',num)) = uicontrol('Style', 'checkbox',...
%         'Units','characters',...
%         'Position', [xpos+10 0 3 1],...
%         'Tag', strcat('showSys',num),...
%         'Parent', handles.ProjectsPanel,...
%         'Callback',eval(strcat('@(hObject,eventdata)MainScreen1(',quote,'popupmenu_Callback',quote,',hObject,eventdata,guidata(hObject))')),...
%         'Visible',vis,'UserData',i);
%         set(handles.(strcat('showSys',num)),'Value',true);
%         %% text show:
%         handles.(strcat('textSys',num)) = uicontrol('Style', 'text', 'String', 'Show:',...
%         'Units','characters',...
%         'Position', [xpos+1 0 8 1],...
%         'Tag', strcat('textSys',num),...
%         'Parent', handles.ProjectsPanel,...
%         'Visible',vis,'UserData',i);
        %% delete button
        handles.(strcat('DeleteSys',num)) = uicontrol('Style', 'pushbutton',...
        'Units','characters',...
        'Position', [xpos+15 0 5 2],...
        'Tag', strcat('DeleteSys',num),....
        'String','X','ForegroundColor','red','BackgroundColor',[.94,.94,.94],'FontSize',18,'FontWeight','bold',...
        'Parent', handles.ProjectsPanel,...
        'Callback',eval(strcat('@(hObject,eventdata)MainScreen1(',quote,'DeleteSys_Callback',quote,',hObject,eventdata,guidata(hObject))')),...
        'Visible',vis,'UserData',i);
    end
end
%make prev and next buttons if they don't exist
if ~isfield(handles,'Copy')
    handles.PrevSys = uicontrol('Style', 'pushbutton', 'String', 'Prev',...
    'Units','characters','FontSize', 10,'Position', [5 0 15 1.8],...
    'Tag', 'PrevSys','Parent', handles.ProjectsPanel,...
    'Callback',eval(strcat('@(hObject,eventdata)MainScreen1(',quote,'PrevSys_Callback',quote,',hObject,eventdata,guidata(hObject))')),...
    'Visible','off','UserData',1);

    handles.NextSys = uicontrol('Style', 'pushbutton', 'String', 'Next',...
    'Units','characters','FontSize', 10,'Position', [70 0 15 1.8],...
    'Tag', 'NextSys','Parent', handles.ProjectsPanel,...
    'Callback',eval(strcat('@(hObject,eventdata)MainScreen1(',quote,'NextSys_Callback',quote,',hObject,eventdata,guidata(hObject))')),...
    'Visible','off','UserData',1);
    if nSys>3 
        set(handles.NextSys,'Visible','on');
    end

    %make copy, load, and save buttons
    handles.Copy = uicontrol('Style', 'pushbutton', 'String', 'Copy',...
    'Units','characters','FontSize', 10,'Position', [85 2.5 15 2],'BackgroundColor',[.76,.87,.78],...
    'Tag', 'Copy','Parent', handles.ProjectsPanel,...
    'Callback',eval(strcat('@(hObject,eventdata)MainScreen1(',quote,'Copy_Callback',quote,',hObject,eventdata,guidata(hObject))')),'Visible','on');

    handles.Load= uicontrol('Style', 'pushbutton', 'String', 'Load',...
    'Units','characters','FontSize', 10,'Position', [85 .5 15 2],'BackgroundColor',[.3,.75,.93],...
    'Tag', 'Load','Parent', handles.ProjectsPanel,...
    'Callback',eval(strcat('@(hObject,eventdata)MainScreen1(',quote,'Load_Callback',quote,',hObject,eventdata,guidata(hObject))')),'Visible','on');

    handles.Save = uicontrol('Style', 'pushbutton', 'String', 'Save',...
    'Units','characters','FontSize', 10,'Position', [4 0 10 2],'BackgroundColor',[.85,.33,.10],...
    'Tag', 'Save','Parent', handles.ProjectsPanel,...
    'Callback',eval(strcat('@(hObject,eventdata)MainScreen1(',quote,'Save_Callback',quote,',hObject,eventdata,guidata(hObject))')),'Visible','on');
end