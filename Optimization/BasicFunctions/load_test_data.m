% This script loads the test data that the projects will be run against
%Option a) use simulated building if anything in Plant.Building (need to double check for other needed data streams like solar/wind/hydro)
%Option b) use historical record already in Plant.Data
%Option c) load test data from a file
global Plant TestData Model_dir
if isempty(TestData)
    if isfield(Plant,'Building') && ~isempty(Plant.Building)
        %use building/s in project
        k = 1;
    elseif isfield(Plant,'Data') && ~isempty(Plant.Data)
        %use data in project
%         k = center_menu('Select Test Data', 'Simulated Building', 'Historical Data in project', 'Load from File');
        k = 3;
    else
        k = center_menu('Select Test Data', 'Simulated Building', 'Load from File');
        if k ==2
            k = 3;
        end
    end
    if k == 1
        if isfield(Plant,'Building') && ~isempty(Plant.Building)
            %use building/s in project
        else
            files = dir(fullfile(Model_dir, 'System Library','Buildings','*.mat'));
            list_b=strrep({files.name},'.mat','');
            [s,v] = listdlg('ListString',list_b,'SelectionMode','multiple','Name','Saved buildings','PromptString','Select one or more buildings');
            Plant.Building = [];
            for i = 1:1:length(s)
                load(fullfile(Model_dir, 'System Library','Buildings',strcat(list_b{s(i)},'*.mat')));
                Plant.Building(i) = building;
            end
        end
        %% select simulation dates
        n_s = 365*24/Plant.optimoptions.Resolution+1;
        TestData.Timestamp = linspace(datenum([2017,1,1]),datenum([2018,1,1]),n_s)';
        %load weather file if not already loaded
        if ~isfield(TestData,'Weather') || isempty(TestData.Weather)
            if ~isfield(Plant,'Weather') || isempty(Plant.Weather)
                cd(fullfile(Model_dir,'System Library','Weather'))
                [fn,pn,~] = uigetfile('*.mat','Load Weather File');
                weather = import_weather(fullfile(pn,fn));
                cd(Model_dir)
                TestData.Weather = interpolate_weather(weather,TestData.Timestamp);
            else
                TestData.Weather = interpolate_weather(Plant.Weather,TestData.Timestamp);
            end
        end
        
    elseif k == 2
        TestData = Plant.Data;%% Revise this so you can pull from more than what is loaded in Plant
        if isfield(Plant.Data,'Hydro')
            TestData = rmfield(TestData,'Hydro');
            TestData.Hydro.SourceSink = Plant.Data.Hydro.SourceSink;
            TestData.Hydro.OutFlow = Plant.Data.Hydro.OutFlow;
        end  
    elseif k == 3
        cd(fullfile(Model_dir,'Data','TestData'))
        [fn,pn,~] = uigetfile('*.mat','Load Test Data File');
        load(fullfile(pn,fn));
        cd(Model_dir)
    end
end
if isfield(Plant,'Building') && ~isempty(Plant.Building) && ~strcmp(Plant.optimoptions.forecast,'Building')
    TestData.Building = load_test_building(Plant.Building,Plant.Network,TestData.Timestamp,TestData.Weather);
end