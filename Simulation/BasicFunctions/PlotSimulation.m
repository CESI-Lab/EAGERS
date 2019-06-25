function PlotSimulation(T,Y,plotFC,animateFC,PlotTurbos)
global TagInf  modelParam LinMod
if ~isempty(modelParam)
    Mod = modelParam;
else Mod = LinMod;
end
time = TagInf.Time;
if isfield(Mod,'Plot')
    n = length(Mod.Plot);
    for i = 1:1:n
        tagName = Mod.Plot{i};
        r = strfind(tagName,'.');
        block = tagName(1:r-1);
        tag = tagName(r+1:end);
        figure(i)
        hold off
        if isfield(TagInf,block) && isfield(TagInf.(block),tag)
            plot(time,TagInf.(block).(tag))
            hold on
            ylabel(tagName);
        end
    end
end
if plotFC
    CompNames = fieldnames(Mod.Components);
    for i = 1:1:length(CompNames)
        if isfield(Mod.Components.(CompNames{i}),'type') && (strcmp(Mod.Components.(CompNames{i}).type,'FuelCell') || strcmp(Mod.Components.(CompNames{i}).type,'Electrolyzer'))
            block = Mod.Components.(CompNames{i});
            CellMap(Y(end,:),block,n+1,'T');
            CellMap(Y(end,:),block,n+2,'I');
            n = n+2;
        end
    end 
end
if animateFC
    CompNames = fieldnames(Mod.Components);
    for i = 1:1:length(CompNames)
        if isfield(Mod.Components.(CompNames{i}),'type') && (strcmp(Mod.Components.(CompNames{i}).type,'FuelCell') || strcmp(Mod.Components.(CompNames{i}).type,'Electrolyzer'))
            block = Mod.Components.(CompNames{i});
            Animate(T,Y,block,n+1)
            n = n+1;
        end
    end
end
if PlotTurbos
    CompNames = fieldnames(Mod.Components);
    for i = 1:1:length(CompNames)
        if strcmp(Mod.Components.(CompNames{i}).type,'Blower') 
            CTmap(Mod.Components.(CompNames{i}),[],n+1)
            n = n+1;
        end
        if strcmp(Mod.Components.(CompNames{i}).type,'Compressor') 
            comp = Mod.Components.(CompNames{i});
            CTmap(comp,[],n+1)
            n = n+1;
        end
        if strcmp(Mod.Components.(CompNames{i}).type,'Turbine') 
            turb = Mod.Components.(CompNames{i});
            CTmap([],turb,n+1)
            n = n+1;
        end
    end
end
end %Ends function PlotSimulation