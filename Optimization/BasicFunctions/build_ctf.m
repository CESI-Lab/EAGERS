function ctf = build_ctf(building)
%Energy balance for surfaces
%dT_surf = A*Tsurf + B*[Tzone;Texterior];
% Q = C*Tsurf + D*[Tzone;Texterior];
%where A = ctf.A, B = ctf.B,...

%absorbed solar + longwave + convection - conduction = 0
%absorbed solar shows up as surface gain
%longwave = h_ground*(Tground-Tsur) + h_sky(T_sky- T_surface) + h_air(T_air-Tsurface); %assume  T_ground = T_air
%convection = h*Area*(T_zone - T_surface)
%conduction = 1/R * (T_boundary - T_surface); %where T_boundary is the other connected surface
%temperature(altitude) = T (weather station) - 0.0065 K/m * altitude
% i_nodes = 5;%interior wall nodes (2 means 1 at each surface, 3+ includes middle ones)
n_zone = length(building.zones.name);
surf = building.surfaces;
n_sur = length(surf.name);
% i_nodes = surf.i_nodes;
ctf.z_height = zeros(n_zone,1);
s_states = 0;%number of temperature states associated with walls
T_ext = 0;
ctf.sur_state1 = ones(n_sur,1);
ctf.sur_state2 = zeros(n_sur,1);
ctf.map_sur2zone = zeros(n_zone,n_sur);
ctf.map_zone2sur = zeros(n_zone,n_sur);
z_names = building.zones.name;
for s = 1:1:n_sur
    j = nonzeros((1:length(z_names))'.*strcmpi(surf.zone{s},z_names));%zone that it is touching
    ctf.map_sur2zone(j,s) = 1;%interior zone that this surface has convective heat transfer with
    if strcmp(surf.type{s},'InternalMass')
        ctf.sur_state1(s) = s_states+1;
        s_states = s_states+1;%adiabatic surface only need 1 state
    else
        ctf.sur_state1(s) = s_states+1;
        i_nodes = length(surf.capacitance{s});
        switch surf.boundary{s}
            case {'Outdoors';'Ground'}
                s_states = s_states+i_nodes;
                T_ext = T_ext + 1;
                ctf.sur_state2(s) = s_states;
            case 'Surface'
                b = nonzeros((1:length(surf.name))'.*strcmpi(surf.object{s},surf.name));
                if b<s
                    s_states = s_states+1;%already counted the intermediate states
                elseif b == s 
                    s_states = s_states+1;%adiabatic surface boundary condition is itself (conduction negated, only need 1 state)
                elseif b>s
                    s_states = s_states+i_nodes-1;
                end
        end
    end
end
n_win = length(building.windows.name);
ctf.map_win2zone = zeros(n_zone,n_win);
ctf.map_zone2win = zeros(n_zone,n_win);
for i = 1:1:n_win
    ctf.map_win2zone(building.windows.zone_index(i),i) = 1;
end


n_ext = nnz(ctf.sur_state2);
for i = 1:1:n_zone
    sur = nonzeros((1:n_sur)'.*(ctf.map_sur2zone(i,:)'));
    win = nonzeros((1:n_win)'.*(ctf.map_win2zone(i,:)'));
    zone_surf_area = sum(surf.area(sur))+sum(building.windows.area(win));
    if ~isempty(sur)
        ctf.map_zone2sur(i,sur) = surf.area(sur)/zone_surf_area;
        height = [];
        for j = 1:1:length(sur)
            if ~isempty(surf.vertices{sur(j)})
                height(end+1) = mean(surf.vertices{sur(j)}(:,3));
            end
        end
    end
    ctf.z_height(i) = median(height);
    if ~isempty(win)
        ctf.map_zone2win(i,win) = building.windows.area(win)/zone_surf_area;
    end
end


ctf.A = zeros(s_states);
ctf.R_pos = zeros(s_states);
ctf.R_neg = zeros(s_states);
ctf.B = zeros(s_states,n_zone+n_sur+1); 
ctf.C = zeros(n_sur,s_states);
ctf.D = zeros(n_sur,n_zone);
ctf.capacitance = zeros(s_states,1);
ctf.interior_surface = zeros(s_states,n_sur);%map surfaces to interior surface T state
ctf.exterior_surface = zeros(s_states,n_ext);%map surfaces to exterior surface T state
ctf.exterior_area = zeros(n_ext,1);
ctf.surf_area = zeros(s_states,1);
ctf.subsurf_names = cell(s_states,1);
ctf.s_height = zeros(n_ext,1);
ctf.ground_cond = zeros(0,1);
i = 0;
s = 0;
while s<n_sur
    s = s+1;%move on to next surface
    j = nonzeros((1:n_zone)'.*ctf.map_sur2zone(:,s));
    if strcmp(surf.type{s},'InternalMass')
        C = sum(surf.capacitance{s});%adiabatic surface (not split into pieces)
        ctf.capacitance(i+1) = C;
        ctf.A(i+1,i+1) = - surf.area(s)/C;%convection to zone air (dT/dt = h*A(Tzone - Tsurface)/C
        ctf.B(i+1,j) = surf.area(s)/C;%convection to zone air (dT/dt = h*A(Tzone - Tsurface)/C
        ctf.C(s,i+1) = -1;%convection to zone air (Q = h*(Tzone - Tsurface)
        ctf.D(s,j) = 1;%convection to zone air (Q = h*(Tzone - Tsurface)
        
        ctf.interior_surface(i+1,s) = 1;
        ctf.subsurf_names(i+1) = surf.name(s);
        ctf.surf_area(i+1) = surf.area(s);
        i = i+1;
    else
        i_nodes = length(surf.capacitance{s});
        R = surf.resistance{s};
        C = surf.capacitance{s};
        if strcmp(surf.boundary{s},'Surface')
            b = nonzeros((1:length(surf.name))'.*strcmpi(surf.object{s},surf.name));
            if b == s %boundary is itself, representing with 1 state
                C = sum(surf.capacitance{s});
            end
        end
        ctf.capacitance(i+1) = C(1);
        
        A = surf.area(s);
        ctf.A(i+1,i+1) = - A/C(1);%convection to zone air (dT/dt = h*A(Tzone - Tsurface)/C
        ctf.B(i+1,j) = A/C(1);%convection to zone air (dT/dt = h*A(Tzone - Tsurface)/C
        
        ctf.C(s,i+1) = -1;%convection to zone air (Q = h*(Tzone - Tsurface)
        ctf.D(s,j) = 1;%convection to zone air (Q = h*(Tzone - Tsurface)
        
        ctf.interior_surface(i+1,s) = 1;
        switch surf.boundary{s}
            case {'Outdoors';'Ground'}
                ctf.R_pos(i+1,i+2) = 1/R(1); %Q_2-->1 = 1/R*(T_2 - T_1)
                ctf.R_neg(i+1,i+1) = -1/R(1);%conduction only from the next node
                for m = 2:i_nodes-1%condution in both directions
                    ctf.R_pos(i+m,i+m-1) = 1/R(m-1);
                    ctf.R_neg(i+m,i+m) = -1/R(m-1) - 1/R(m);
                    ctf.R_pos(i+m,i+m+1) = 1/R(m);
                end
                if strcmp(surf.boundary{s},'Ground')
                    ctf.R_pos(i+i_nodes,i+i_nodes-1) = 1/R(i_nodes-2);
                    ctf.R_neg(i+i_nodes,i+i_nodes) = -1/R(i_nodes-2);% -2/R: exterior temperature is appended to list of zone temperature states and multiplied by B, the last temperature state associated with this surface(the side facing the exterior)is part of Tsurf
                    ctf.ground_cond(end+1,1) = 1/R(i_nodes-1);
%                     ctf.B(i+i_nodes,n_zone+s) = 1/(R*C);%conduction to ground
                elseif strcmp(surf.boundary{s},'Outdoors')
                    ctf.R_pos(i+i_nodes,i+i_nodes-1) = 1/R(i_nodes-1);
                    ctf.R_neg(i+i_nodes,i+i_nodes) = -1/R(i_nodes-1);%conduction only from the previous node
                    
                    ctf.A(i+i_nodes,i+i_nodes) = - A/C(i_nodes);%exterior temperature is appended to list of zone temperature states and multiplied by B, the last temperature state associated with this surface(the side facing the exterior)is part of Tsurf
                    ctf.B(i+i_nodes,n_zone+s) = A/C(i_nodes);%convection to air outside zone
                    ctf.B(i+i_nodes,end) = A/C(i_nodes);%convection to sky temperature
                    ctf.ground_cond(end+1,1) = 0;
                end
                ctf.capacitance(i+1:i+i_nodes) = C;
                ext = nnz(ctf.sur_state2(1:s));
                ctf.exterior_surface(i+i_nodes,ext) = 1;
                ctf.exterior_area(ext) = building.surfaces.area(s);
                for s2 = 1:i_nodes
                    ctf.subsurf_names(i+s2) = {strcat(surf.name{s},'_',num2str(s2))};
                end
                ctf.s_height(ext) = mean(surf.vertices{s}(:,3));
                ctf.surf_area(i+1:i+i_nodes) = surf.area(s);
                i = i+i_nodes;
            case 'Surface'
                if b>s%this avoids doing this twice, since calculations for surface A-->B and B-->A are the same
                    ctf.R_pos(i+1,i+2) = 1/R(1);
                    ctf.R_neg(i+1,i+1) = - 1/R(1);
                    for m = 2:i_nodes-2%condution in both directions
                        ctf.R_pos(i+m,i+m-1) = 1/R(m-1);
                        ctf.R_neg(i+m,i+m) = -1/R(m-1) - 1/R(m);
                        ctf.R_pos(i+m,i+m+1) = 1/R(m);
                    end
                    m = i_nodes-1;
                    ctf.R_pos(i+m,i+m-1) = 1/R(m-1);
                    ctf.R_neg(i+m,i+m) = -1/R(m-1) - 1/R(m);
                    ctf.R_pos(i+m,ctf.sur_state1(b)) = 1/R(m);
                    
                    ctf.R_pos(ctf.sur_state1(b),i+m) = 1/R(i_nodes-1);
                    ctf.R_neg(ctf.sur_state1(b),ctf.sur_state1(b)) = -1/R(i_nodes-1);
                    ctf.capacitance(i+1:i+i_nodes-1) = C(1:i_nodes-1);
                    for s2 = 1:i_nodes-1
                        ctf.subsurf_names(i+s2) = {strcat(surf.name{s},'_',num2str(s2))};
                    end
                    ctf.surf_area(i+1:i+i_nodes-1) = surf.area(s);
                    i = i+i_nodes-1;
                else
                    %adiabatic surface when b == s, surface boundary condition is itself (conduction negated)
                    %when b<s, already did all intermediate temperatures for this wall
                    ctf.subsurf_names(i+1) = surf.name(s);
                    ctf.surf_area(i+1) = surf.area(s);
                    i = i+1;
                end
        end
    end
end
ctf.sur_state2 = nonzeros(ctf.sur_state2);
end%Ends function build_ctf