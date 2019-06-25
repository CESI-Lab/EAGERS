function plotLines(GenDisp,Time)
%%plot line power transfers and losses
global Plant
networkNames = fieldnames(Plant.Network);
networkNames = networkNames(~strcmp('name',networkNames));
networkNames = networkNames(~strcmp('Equipment',networkNames));
networkNames = networkNames(~strcmp('Location',networkNames));
networkNames = networkNames(~strcmp('Buildings',networkNames));

dt = Time' - [0, Time(1:end-1)]';
horizon = Plant.optimoptions.Horizon;
if horizon>24
    axStep = floor(horizon/10);
elseif horizon>12
    axStep=2;
else
    axStep=1;
end
axTick = ceil(Time(1)-dt(1)):axStep:Time(end);
axIndex = axTick;
while axIndex(end)>24
    axIndex(axIndex>24) = axIndex(axIndex>24)-24;
end

LN = Plant.subNet.lineNames;
n = length(LN);
LV = GenDisp(2:end,nG+1:end);%Power Transfer across all lines
LLN = zeros(nS,nL);
for t = 1:1:nS
    LLN = LV(t,:).*LineLoss;
end
for q = 1:1:length(networkNames)
    names = {};
    Transfer =[];
    for i = 1:1:n
        name = LN{i};
        a = strfind(name,'_');
        net = name(a(1)+1:a(2)-1);
        if strcmp(net,networkNames{q})
            Transfer(1:nS,end+1) = LV(:,i);
        end 
        names(end+1) = {strcat(name(1:a(1)),name(a(2)+1:end))};
    end 
    %% Plot
    if ~isempty(Transfer)
        figure(q+1)
        hold off
        plotTime = sort([Time(1)-1,Time,Time-.9999]);%Time(end)+dt(end)-.0001]);%make sure that there is a step wise characteristic of the plot
        h1 = plot(plotTime,Transfer);
        if strcmp(networkNames{q},'DistrictCool')
            colormap cool
        elseif strcmp(networkNames{q},'DistrictHeat')
            colormap autumn
        else
            colormap colorcube
            LTColors = {[.6 0 .6],[1 .8 .2],[1 0 0],[.2 1 1],[.4 0 .6]};%purple [.4 0 .6] if using jet, [.6 0 .6] if default map, orange, cyan
            LLColors = {[.9 0 .9],[1 .9 .2],[1 .2 0],[.5 1 1],[.7 0 .9]};
            for i = 1:1:n
                if i > length(LN)
                    set(h1(i),'FaceColor',LLColors{i-length(LN)});
                else
                    set(h1(i),'FaceColor',LTColors{i});
                end 
            end
        end
        legend(names,'Fontsize',16)
%         hold on
%         if ~isempty(negBars)
%             OoM = log10(max(sum(Transfer,2)+sum(negBars,2)));
%         else OoM = log10(max(sum(Transfer,2)));
%         end
%         if (OoM-floor(OoM))==0 %count in increments of 1, 10, 100 or 1000 etc
%             Yspace = 10^(OoM-1);
%             Ymax = 10^OoM;
%         elseif (OoM-floor(OoM))> 0.6990 %count in increments of 1, 10, 100 or 1000 etc
%             Yspace = 10^floor(OoM);
%             Ymax = 10^ceil(OoM);
%         elseif (OoM-floor(OoM))> 0.30103 %count in increments of 5, 50, 500 or 5000 etc
%             Yspace = .5*10^floor(OoM);
%             Ymax = .5*10^ceil(OoM);
%         else  %count in increments of 2, 20, 200 or 2000 etc
%             Yspace = .2*10^floor(OoM);
%             Ymax = .2*10^ceil(OoM);
%         end
%         negTicks = 0;
%         Ymin = Yspace*negTicks;
%         if isempty(Ymin)
%             Ymin = 0;
%         end
%         pTicks = Ymax/Yspace;
%         legend(name,'Fontsize',16,'Orientation','Horizontal','Location','NorthOutside','Box','off')
%         legend([name name2],'FontSize',16)
        xlim([Time(1)-dt(1), Time(end)])%[0,25]
%         ylim([Ymin,Ymax])
%         set(gca,'YTick',Ymin:Yspace:Ymax,'FontSize',14)
%         set(gca, 'XTick',axTick,'XTickLabel', {axIndex})%{xticks(1:nS+1)})
        xlabel('Time (hour)','Color','k','FontSize',18)
        ylabel('Power','Color','k','FontSize',18)
    end
end 