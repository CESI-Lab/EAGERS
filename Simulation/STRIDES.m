%% STRIDES %%
% This script gives the option to initialize or load a saved non-linear systems model
% Next it offers to linearize that model,or load a saved linearization. 
% Finally it is able to run a transient on either/both the linear or non-linear model
global Model_dir SimSettings WaitBar modelParam LinMod IterCount TagInf TagFinal Tags Inlet Outlet
%%clear variables
LinMod =[];
modelParam = [];
TagInf = [];
TagFinal = [];
Tags = [];
SimSettings = [];
Inlet = [];
Outlet = [];
%% Load a Plant
%%%%-- either load a plant and initialize it, or load a saved Plant & model Parameters
%%%%-- Working options are: SOFCstack, SOECstack, SOFCsystem, GasTurbine
J = center_menu('Non-linear model options','Initialize from system description','Load previously initialized plant','Skip to pre-loaded linear model');
if J ==1
    %the goal is to have a GUI replace these m-files
    ModelFiles=dir(fullfile(Model_dir,'Model_Library','*.m'));
    list = {};
    for i = 1:1:length(ModelFiles)
        list(i) = cellstr(strrep(ModelFiles(i).name,'.m',''));
    end
    [s,OK] = listdlg('PromptString','Select Model', 'SelectionMode','single','ListString',list);
    if OK~=1
        disp('Invalid selection. Exiting...')
    else
        modelName = list{s};
        Plant = feval(modelName);
        WaitBar.Show = 1;
        BuildModel(Plant); %%Build & Initialize model
        J2 = center_menu('Save Model?','Yes','No');
        if J2 ==1
            [f,p]=uiputfile(fullfile(Model_dir,'Model_Library','Saved_Models',strcat(modelName,'.mat')),'Save Model As...');
            save([p f],'modelParam')
        end
    end
elseif J ==2
    ModelFiles=dir(fullfile(Model_dir,'Model_Library','Saved_Models','*.mat'));
    list = {};
    for i = 1:1:length(ModelFiles)
        list(i) = cellstr(ModelFiles(i).name);
    end
    [s,OK] = listdlg('PromptString','Select Model', 'SelectionMode','single','ListString',list);
    if OK~=1
        disp('Invalid selection. Exiting...')
    else
        modelName = list{s};
        modelName = strrep(modelName,'.mat','');
        load(fullfile(Model_dir,'Model_Library','Saved Models',modelName));
        Outlet = modelParam.NominalOutlet; SimSettings = modelParam.NominalSettings; Tags = modelParam.NominalTags;
    end
end
%% Create or load a set of linear models
if J ==3
    J2 = 2; %load linearization
else J2 = center_menu('Linear model option','Create new linearization','Load previously linearized plant','Skip linearization');
end
if J2 ==1
    CreateLinModel;
    J3 = center_menu('Save Linearized Model?','Yes','No');
    if J3 ==1
        [f,p]=uiputfile(fullfile(Model_dir,'Model_Library','Saved Linearizations',strcat(modelName,'.mat')),'Save Linearized Model As...');
        save([p f],'LinMod')
    end
elseif J2 ==2
    ModelFiles=dir(fullfile(Model_dir,'Model_Library','Saved Linearizations','*.mat'));
    list = {};
    for i = 1:1:length(ModelFiles)
        list(i) = cellstr(ModelFiles(i).name);
    end
    [s,OK] = listdlg('PromptString','Select Model', 'SelectionMode','single','ListString',list);
    if OK~=1
        disp('Invalid selection. Exiting...')
    else
        load(fullfile(Model_dir,'Model_Library','Saved Linearizations',list{s}));
    end
end

if J2 ==3
    if J ==1 || J ==2
        J3 = center_menu('Simulation Options','Simulate non-linear response','Exit');
        if J3 == 2
            J3 = 4;%exit
        end
    else
        J3 = 4; %no linear or non-linear model loaded, exit
    end
elseif J == 3
    J3 = 2; %simulate linear response
else
    J3 = center_menu('Simulation Options','Simulate non-linear response','Simulate linear response','Simulate both linear and non-linear response','Neither');
end
if J3 ~=4 
    %set up the transient
    %first identify any controller input (lookup functions), let user pick ones with a schedule
    %then get the variables from that function and allow the user to edit them
    if ~isempty(modelParam)
        controls = fieldnames(modelParam.Controls);
        Outlet = modelParam.NominalOutlet; SimSettings = modelParam.NominalSettings; Tags = modelParam.NominalTags;
    elseif ~isempty(LinMod)
        controls = fieldnames(LinMod.Controls);
        Outlet = LinMod.NominalOutlet; SimSettings = LinMod.NominalSettings; Tags = LinMod.NominalTags;
    end
    A = (inputdlg({'Test Duration (hr)';'Maximum Step Size (hr)'},'Specify length of the transient simulation in hours',1,{'24';'.25';}));
    SimSettings.RunTime = eval(A{1})*3600;
    options = odeset('MaxStep',eval(A{2})*3600);
    t_Steps = [0,SimSettings.RunTime];
    for i = 1:1:length(controls)
        if ~isempty(modelParam)
            Cont = modelParam.Controls.(controls{i});
        elseif ~isempty(LinMod)
            Cont = LinMod.Controls.(controls{i});
        end
        for j = 1:1:length(Cont.connections)
            r = strfind(Cont.connections{j},'.');
            if ~isempty(Cont.connections{j}) && ~isnumeric(Cont.connections{j}) && isempty(r) %must be a lookup function
                [globvar,Prompt,DefaultVal] = feval(Cont.connections{j},0,'loadparam'); 
                A = inputdlg(Prompt,strcat('Specify the transient input parameters for the function',Cont.connections{j},'Any string will be evaluated, but must create vectors of equal length.'),1,DefaultVal);
                for k = 1:1:length(globvar)
                    SimSettings.(globvar{k}) = eval(A{k});
                end
                t_Steps = [t_Steps SimSettings.(globvar{1})*3600];
            end
        end
        SimSettings.t_Steps = sort(unique(t_Steps));
        %% modify control terms
        nC = length(Cont.Gain);
        Prompt = {};
        DefaultVal = {};
        for j = 1:1:nC
            Prompt(end+1) = cellstr(strcat(Cont.PIdescription{j},'---Proportional'));
            Prompt(end+1) = cellstr(strcat(Cont.PIdescription{j},'---Integral'));
            DefaultVal(end+1) = cellstr(num2str(Cont.PropGain(j)));
            DefaultVal(end+1) = cellstr(num2str(Cont.Gain(j)));
        end
        A= inputdlg(Prompt,'Specify the controller parameters',1,DefaultVal);
        for j = 1:1:nC
            Cont.PropGain(j) = str2double(A(2*j-1));
            Cont.Gain(j) = str2double(A(2*j));
        end
        if ~isempty(modelParam)
            modelParam.Controls.(controls{i}) = Cont;
        end
        if ~isempty(LinMod)
            LinMod.Controls.(controls{i}) = Cont;
        end
    end 
end

%% Run a transient on non-linear model
if J3 ==1 || J3 == 3
    WaitBar.Show = 1; 
    IterCount = 1; TagInf =[]; TagFinal =[];  WaitBar.Text = 'Running non-linear model with transient';WaitBar.Handle =waitbar(0,WaitBar.Text);
    tic; [T, Y] = ode15s(@RunBlocks, [0, SimSettings.RunTime], modelParam.IC,options); disp(strcat('Time to run model:',num2str(toc),' seconds'));close(WaitBar.Handle);
    PlotSimulation(T,Y,1,0,1)% Plot the tags and scopes in Plant.Plot. The first option (after Y) is to plot any fuel cell or electrolyzer temperature profiles, the second option is to mave a video of the transient, the third is to plot compressor turbine and blower maps
end

%% Run transient on set of linear models
if J3 ==2 || J3 ==3
    % need to find better initial condition (won't always start at nominal power)
    IC =  [LinMod.Model{1}.X0;LinMod.Model{1}.UX0];
    IterCount = 1; TagInf =[]; TagFinal =[]; WaitBar.Show = 1; WaitBar.Text = 'Running linear model with transient';WaitBar.Handle =waitbar(0,WaitBar.Text);
    tic; [T, Y] = ode15s(@RunLinSystem, [0, SimSettings.RunTime], IC,options); disp(strcat('Time to run model:',num2str(toc),' seconds'));close(WaitBar.Handle);
    if J3 == 2
        PlotSimulation(T,Y,1,0,1)
    else
%         PlotSimulation(T,Y,0,0,0) % have already plotted the maps, just adding linear response to non-linear response
    end
end