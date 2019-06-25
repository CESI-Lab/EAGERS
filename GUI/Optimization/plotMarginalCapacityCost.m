function plotMarginalCapacityCost(handles,market_margin_cost,date)
global Plant
axes(handles.MarketAxes)
% cla(handles.MarketAxes)
list = get(handles.MarketPopup,'String');
str =  list{get(handles.MarketPopup,'Value')};
set(handles.CapacityTable,'Data',Plant.Market.Price);
switch str
    case 'Market Prices'
        set(handles.MarketTitle,'String','Capacity Market Prices ($/kW)');
        [X,I] = sort([linspace(0,23,24), linspace(.9999,22.9999,23),24]);
        Y = [Plant.Market.Price(1,:),Plant.Market.Price(1,:)];
        Y = Y(I);
        hold(handles.MarketAxes,'off')
        plot(handles.MarketAxes,X,Y,'k')
        hold(handles.MarketAxes,'on')
        Y = [Plant.Market.Price(2,:),Plant.Market.Price(2,:)];
        Y = Y(I);
        plot(handles.MarketAxes,X,Y,'b')
        Y = [Plant.Market.Price(3,:),Plant.Market.Price(3,:)];
        Y = Y(I);
        plot(handles.MarketAxes,X,Y,'r')
        ylabel('Price ($/kW)')
        legend('Day Ahead','Hour Ahead','Realtime')
    case 'Reserve Capacity'
        set(handles.MarketTitle,'String','On-Site Capacity (kW)');
        D = floor(date)+linspace(0,1,24);
        D(D<market_margin_cost.Timestamp(1)) = D(D<market_margin_cost.Timestamp(1))+1;
        D(D>market_margin_cost.Timestamp(end)) = market_margin_cost.Timestamp(end);
        [X,I] = sort([linspace(0,23,24), linspace(.9999,22.9999,23),24]);
        Y = interp1(market_margin_cost.Timestamp,Plant.Market.Capacity(1,:),D); Y = [Y,Y];Y = Y(I);
        hold(handles.MarketAxes,'off')
        plot(handles.MarketAxes,X,Y,'k')
        hold(handles.MarketAxes,'on')
        Y = interp1(market_margin_cost.Timestamp,Plant.Market.Capacity(2,:),D); Y = [Y,Y];Y = Y(I);
        plot(handles.MarketAxes,X,Y,'b')
        Y = interp1(market_margin_cost.Timestamp,Plant.Market.Capacity(3,:),D); Y = [Y,Y];Y = Y(I);
        plot(handles.MarketAxes,X,Y,'r')
        ylabel('Capacity (kW)')
        legend('Generator Capacity','Demand Response Capacity','Committed Capacity')
    case {'Generation ($/kW)';'Demand Response ($/kW)'}
        if strcmp(str,'Generation ($/kW)')
            set(handles.MarketTitle,'String','Marginal Capacity Cost ($/kW)');
            f1 = 'GenCapacity';
            f2 = 'GenCost';
        else
            set(handles.MarketTitle,'String','Marginal DR Cost ($/kW)');
            f1 = 'DR_Capacity';
            f2 = 'DR_Cost';
        end
        step = 10;
        n = floor(max(max(market_margin_cost.(f1)/step)));
        Z = zeros(n+1,25);
        D = floor(date) + linspace(1,24,24)/24;
        D(D<market_margin_cost.Timestamp(1)) = D(D<market_margin_cost.Timestamp(1))+1;
        D(D>market_margin_cost.Timestamp(end)) = market_margin_cost.Timestamp(end);
        for i = 1:1:24
            k = nnz(market_margin_cost.Timestamp<=D(i));
            if k == length(market_margin_cost.Timestamp)
                k = k-1;
            end
            r = (D(i) - market_margin_cost.Timestamp(k))/(market_margin_cost.Timestamp(k+1) - market_margin_cost.Timestamp(k)); %weighting of next interpolation point
            n2 = floor(min(max(market_margin_cost.(f1)(:,k)/step),max(market_margin_cost.(f1)(:,k+1)/step)));%max kW of spare gen capacity
            if n2>0
                Y = market_margin_cost.(f2)(:,k)./market_margin_cost.(f1)(:,k);
                Z1 = interp1([0;market_margin_cost.(f1)(:,k)],[Y(1);Y],linspace(0,n2*step,n2+1)');
                Y = market_margin_cost.(f2)(:,k+1)./market_margin_cost.(f1)(:,k+1);
                Z2 = interp1([0;market_margin_cost.(f1)(:,k+1)],[Y(1);Y],linspace(0,n2*step,n2+1)');
                Z(1:n2+1,i) = (1-r)*Z1+r*Z2;
                Z(n2+2:end,i) = NaN;
            else
                Z(:,i) = NaN;
            end            
        end
        Z(:,25) = NaN;
        hold(handles.MarketAxes,'off')
        if any(any(~isnan(Z)))
            pcolor(handles.MarketAxes,[0:1:24],[0:1:n],Z)
            c = colorbar('peer',handles.MarketAxes);
        end
        if strcmp(str,'Generation ($/kW)')
            ylabel('Excess Generation Capacity (kW)')
        else
            ylabel('Demand Response Capacity (kW)')
        end
        % c.Label.String = 'Marginal Cost ($/kW)'; 
end
xlabel('Hour of day')
xlim([0,24])
axTick = 0:2:24;
axIndex = mod(axTick,24);
axIndex([false,axIndex(2:end)==0]) = 24;
set(handles.MarketAxes,'XTick',axTick,'XTickLabel', {axIndex})
hold(handles.MarketAxes,'off')
end %Ends plotMarginalCapacityCost