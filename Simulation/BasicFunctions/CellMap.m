function CellMap(Y,block,FigNum,mode)
if ~isfield(block,'Manifold')
    block.Manifold = 0;
end
columns = block.columns;
rows = block.rows;
nodes = block.nodes;
switch block.direction
    case 'coflow'
        Xlabel = 'Fuel/Air Flow Direction -->';
        Ylabel = 'Cell Width (parallel channels)';
    case 'counterflow'
        Xlabel = 'Fuel Flow --> & Air Flow <--';
        Ylabel = 'Cell Width (parallel channels)';
    case 'crossflow'
        Xlabel = 'Cell Length (m)';
        Ylabel = 'Cell Width (m)';
end

wC = block.L_Cell/columns;
wR = block.W_Cell/rows;
FCstates = Y(block.States).*block.Scale';
if strcmp(mode,'T')
    CellProperty = FCstates(2*nodes+1:3*nodes)-273;%electrolyte temperature
else
    CellProperty = FCstates(end-nodes-1:end-2)/(wC*wR*1e4);%current density
end
if rows ==1
    X = [0 linspace(.5*wC,block.L_Cell -.5*wC,columns) block.L_Cell];
    Y = [0 block.W_Cell];
    Z(1,1:columns+2) = 0;
    C(1,1:columns+2) =  [CellProperty(1) CellProperty(1:columns) CellProperty(columns)];
    Z(2,1:columns+2) = 0;
    C(2,1:columns+2) =  [CellProperty(1) CellProperty(1:columns) CellProperty(columns)];
    %Draw line coordinates
    for j = 1:1:rows+1
        HorizLinesX(3*j-2:3*j) = [0 block.L_Cell 0];
        HorizLinesY(3*j-2:3*j) = (j-1)*wR;
        HorizLinesZ(3*j-2:3*j) = 0;
    end
    for j = 1:1:columns+1
        VertLinesX(3*j-2:3*j) =(j-1)*wC;
        VertLinesY(3*j-2:3*j) = [0 block.W_Cell 0];
        VertLinesZ(3*j-2:3*j) = 0;
    end
    %draw box coordinates
    Box1X = [0 block.L_Cell block.L_Cell 0 0];
    Box1Y = [0 0 block.W_Cell block.W_Cell 0];
    Box1Z = [0 0 0 0 0];
elseif block.Manifold== 0 && rows>1
    %Organize Temperature & node coordinates
    X = [0 linspace(.5*wC,block.L_Cell -.5*wC,columns) block.L_Cell];
    Y = [0 linspace(.5*wR,block.L_Cell -.5*wR,rows) block.W_Cell];
    for j = 1:1:rows
        Z(j+1,2:columns+1) = 0;
        C(j+1,2:columns+1) = CellProperty(1+(j-1)*columns:j*columns);
    end
    Z(1,2:columns+2) = 0;
    C(1,2:columns+2) =  [CellProperty(1:columns) CellProperty(columns)];
    Z(2:rows+2,columns+2) = 0;
    C(2:rows+2,columns+2) = [CellProperty(linspace(columns,rows*columns,rows)) CellProperty(rows*columns)];
    Z(rows+2,linspace(columns+1,1,columns+1)) = 0;
    C(rows+2,linspace(columns+1,1,columns+1)) = [CellProperty(linspace(rows*columns,(rows-1)*columns+1,columns)) CellProperty((rows-1)*columns+1)];
    Z(linspace(rows+1,1,rows+1),1) = 0;
    C(linspace(rows+1,1,rows+1),1) = [CellProperty(linspace((rows-1)*columns+1,1,rows)) CellProperty(1)];
    %Draw line coordinates
    for j = 1:1:rows+1
        HorizLinesX(3*j-2:3*j) = [0 block.L_Cell 0];
        HorizLinesY(3*j-2:3*j) = (j-1)*wR;
        HorizLinesZ(3*j-2:3*j) = 0;
    end
    for j = 1:1:columns+1
        VertLinesX(3*j-2:3*j) =(j-1)*wC;
        VertLinesY(3*j-2:3*j) = [0 block.W_Cell 0];
        VertLinesZ(3*j-2:3*j) = 0;
    end
    %draw box coordinates
    Box1X = [0 block.L_Cell block.L_Cell 0 0];
    Box1Y = [0 0 block.W_Cell block.W_Cell 0];
    Box1Z = [0 0 0 0 0];
elseif block.Manifold == 1  %% need to edit this
    EdgeT = EdgeStates(length(Time),(Slice-1)*2*(rows+columns+2)+1:Slice*2*(rows+columns+2));
    wE1 = (block.L_Plate-block.L_Cell)/2;
    wE2 = (block.W_Plate-block.W_Cell)/2;
    %Organize Temperature & node coordinates
    X = [0 .5*wE1 linspace(wE1+.5*wC,block.L_Plate -(wE1+.5*wC),columns) block.L_Plate-.5*wE1 block.L_Plate];
    Y = [0 .5*wE2 linspace(wE2+.5*wR,block.W_Plate -(wE2+.5*wR),rows) block.W_Plate-.5*wE2 block.W_Plate];
    for j = 1:1:rows
        Z(j+2,3:columns+2) = 0;
        C(j+2,3:columns+2) = CellProperty(1+(j-1)*columns:j*columns);
    end
    %Edges
    Z(2,3:columns+3) = 0;
    C(2,3:columns+3) =  EdgeT(1:columns+1);
    Z(3:rows+3,columns+3) = 0;
    C(3:rows+3,columns+3) = EdgeT(columns+2:rows+columns+2);
    Z(rows+3,linspace(columns+2,2,columns+1)) = 0;
    C(rows+3,linspace(columns+2,2,columns+1)) = EdgeT(rows+columns+3:rows+2*columns+3);
    Z(linspace(rows+2,2,rows+1),2) = 0;
    C(linspace(rows+2,2,rows+1),2) = EdgeT(rows+2*columns+4:2*rows+2*columns+4);

    Z(1,2:columns+4) = 0;
    C(1,2:columns+4) =  [EdgeT(1) EdgeT(1:columns+1) EdgeT(columns+1)];
    Z(2:rows+4,columns+4) = 0;
    C(2:rows+4,columns+4) = [ EdgeT(columns+2) EdgeT(columns+2:rows+columns+2) EdgeT(rows+columns+2)];
    Z(rows+4,linspace(columns+3,1,columns+3)) = 0;
    C(rows+4,linspace(columns+3,1,columns+3)) = [EdgeT(rows+columns+3) EdgeT(rows+columns+3:rows+2*columns+3) EdgeT(rows+2*columns+3)];
    Z(linspace(rows+3,1,rows+3),1) = 0;
    C(linspace(rows+3,1,rows+3),1) = [EdgeT(rows+2*columns+4) EdgeT(rows+2*columns+4:2*rows+2*columns+4) EdgeT(2*rows+2*columns+4)];

    %Draw line coordinates
    for j = 1:1:rows+3
        HorizLinesX(3*j-2:3*j) = [0 block.L_Plate 0];
        HorizLinesY(3*j-2:3*j) = wE2+(j-2)*wR;
        HorizLinesZ(3*j-2:3*j) = 0;
    end
    HorizLinesY(1:3) = 0;
    HorizLinesY(3*j-2:3*j) = block.W_Plate;
    for j = 1:1:columns+3
        VertLinesX(3*j-2:3*j) = wE1+(j-2)*wC;
        VertLinesY(3*j-2:3*j) = [0 block.W_Plate 0];
        VertLinesZ(3*j-2:3*j) = 0;
    end
    VertLinesX(1:3) = 0;
    VertLinesX(3*j-2:3*j) = block.L_Plate;
    %draw box coordinates
    Box1X = [0 block.L_Plate block.L_Plate 0 0];
    Box1Y = [0 0 block.W_Plate block.W_Plate 0];
    Box1Z = [0 0 0 0 0];
    Box2X = [wE1 wE1+block.L_Cell wE1+block.L_Cell wE1 wE1];
    Box2Y = [wE2 wE2 wE2+block.W_Cell wE2+block.W_Cell wE2];
    Box2Z = [0 0 0 0 0];
end

figure(FigNum);
surf(X, Y, Z, C,'FaceColor','interp','EdgeColor','none')
if strcmp(mode,'T')
    caxis([block.TpenAvg-.75*block.deltaTStack-273 block.TpenAvg+.75*block.deltaTStack-273]);
else
    caxis([0 1.5]);%CellProperty
end
hold on
plot3(HorizLinesX,HorizLinesY,HorizLinesZ,'k-','LineWidth',1);
plot3(VertLinesX,VertLinesY,VertLinesZ,'k-','LineWidth',1)
plot3(Box1X,Box1Y,Box1Z,'k-','LineWidth',4)
if block.Manifold ==1
    plot3(Box2X,Box2Y,Box2Z,'k-','LineWidth',3)
end
c = colorbar;
hold off
if block.Manifold ==0
    axis([0 block.L_Cell 0 block.W_Cell])
else
    axis([0 block.L_Plate 0 block.W_Plate])
end
view([0 90])
xlabel(Xlabel,'FontSize',14)
ylabel(Ylabel,'FontSize',14)
% c.Label.String = 'Temperature (C)';
end%Ends function CellMap