function plotGenEfficiency(Gen,handles)
axes(handles.EffCurve)
hold off
str = {};
if isfield(Gen.Output,'Electricity') && nnz(Gen.Output.Electricity) > 0
    c = Gen.Output.Capacity./Gen.Output.Electricity;
    [AX, H1, H2] = plotyy(Gen.Output.Capacity,Gen.Output.Electricity, ...
        Gen.Output.Capacity(2:end),c(2:end));
    hold on
    set(H1,'Color','k','LineStyle','-','LineWidth',2,'Marker','o')
    str(end+1) = {'Electric'};
end
if isfield(Gen.Output,'DirectCurrent') && nnz(Gen.Output.DirectCurrent) > 0
    c = Gen.Output.Capacity./Gen.Output.DirectCurrent;
    [AX, H1, H2] = plotyy(Gen.Output.Capacity,Gen.Output.DirectCurrent, ...
        Gen.Output.Capacity(2:end),c(2:end));
    hold on
    set(H1,'Color','k','LineStyle','-','LineWidth',2,'Marker','o')
    str(end+1) = {'DirectCurrent'};
end
if isfield(Gen.Output,'Heat') && nnz(Gen.Output.Heat) > 0
    if strcmp(Gen.Type,'CHP Generator')
       plot(AX(1),Gen.Output.Capacity,Gen.Output.Heat,'r-o');
    else
        c = Gen.Output.Capacity./Gen.Output.Heat;
        [AX, H1, H2] = plotyy(Gen.Output.Capacity,Gen.Output.Heat, ...
            Gen.Output.Capacity(2:end),c(2:end));
        set(H1,'Color','r','LineStyle','-','LineWidth',2,'Marker','o')
        hold on
    end
    str(end+1) = {'Heat'};
end
if isfield(Gen.Output,'Steam') && nnz(Gen.Output.Steam) > 0
    c = Gen.Output.Capacity./Gen.Output.Steam; 
    [AX, H1, H2] = plotyy(Gen.Output.Capacity,Gen.Output.Steam, ...
        Gen.Output.Capacity(2:end),c(2:end));
    set(H1,'Color','m','LineStyle','-','LineWidth',2,'Marker','o')
    hold on
    str(end+1) = {'Steam'};
end
if isfield(Gen.Output,'Cooling') && nnz(Gen.Output.Cooling) > 0
    c = Gen.Output.Capacity./Gen.Output.Cooling; 
    [AX, H1, H2] = plotyy(Gen.Output.Capacity,Gen.Output.Cooling, ...
        Gen.Output.Capacity(2:end),c(2:end));
    set(H1,'Color','b','LineStyle','-','LineWidth',2,'Marker','o')
    str(end+1) = {'Cooling'};
end
set(AX,{'ycolor'},{'k';'k'})
set(H2,'Color','g','LineStyle',':','LineWidth',3)
str(end+1) = {'Cost Curve'};

if isfield(Gen.Output,'Cooling') && max(Gen.Output.Cooling) > 1
    ylim([0 max(Gen.Output.Cooling)])
    a = round(max(Gen.Output.Cooling))+1;
    ylim(AX(1),[0,a])
    ylim(AX(2),[0,1])
else
    a = round(max(c))+1;
    ylim(AX(1),[0,1])
    ylim(AX(2),[0,a])
    set(AX(1),'YTick',0:.1:1)
    set(AX(2),'YTick',0:a/10:a)
end
set(get(AX(1),'Ylabel'),'String','Efficiency')
set(get(AX(2),'Ylabel'),'String','Cost Curve Shape')
xlabel(AX(1),'% of Capacity')
legend(str);
title('Efficiency / Cost')
set(handles.EffCurve,'UserData',AX)
set(handles.EffCurve,'Tag','EffCurve')
set(AX(2),'Tag','EffCurveAx2')