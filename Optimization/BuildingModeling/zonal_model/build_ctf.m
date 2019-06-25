function ctf = build_ctf(z_names,surfaces,windows)
%% organize surface and interior surface states for heat transfer calculations
n_zone = length(z_names);
n_sur = length(surfaces.name);
ctf.z_height = zeros(n_zone,1);
s_states = 0;%number of temperature states associated with walls
ctf.sur_state1 = ones(n_sur,1);
ctf.sur_state2 = zeros(n_sur,1);
ctf.map_sur2zone = zeros(n_zone,n_sur);
ctf.map_zone2sur = zeros(n_zone,n_sur);

n_win = length(windows.name);
ctf.map_win2zone = zeros(n_zone,n_win);
ctf.map_zone2win = zeros(n_zone,n_win);
for i = 1:1:n_win
    ctf.map_win2zone(windows.zone_index(i),i) = 1;
end

for i = 1:1:n_zone
    sur = f_index(i,surfaces.zone_index);
    win = f_index(i,windows.zone_index);
    zone_surf_area = sum(surfaces.area(sur))+sum(windows.area(win));
    ctf.map_zone2win(i,win) = windows.area(win)/zone_surf_area;
    ctf.map_zone2sur(i,sur) = surfaces.area(sur)/zone_surf_area;
    height = [];
    for j = 1:1:length(sur)
        if ~isempty(surfaces.vertices{sur(j)})
            height(end+1) = mean(surfaces.vertices{sur(j)}(:,3));
        end
    end
    ctf.z_height(i) = median(height);
end
first_side = true(n_sur,1);
for s = 1:1:n_sur
    switch surfaces.boundary{s}
        case {'Outdoors';'Ground'}
            s_states = s_states + length(surfaces.capacitance{s});
        case 'Surface'
            b = nonzeros((1:length(surfaces.name))'.*strcmpi(surfaces.object{s},surfaces.name));
            if b<s
                first_side(s) = false;%already counted the states, this avoids doing things twice, since calculations for surface A-->B and B-->A are the same
            else
                s_states = s_states + length(surfaces.capacitance{s});
            end
    end
end

ctf.capacitance = zeros(s_states,1);
ctf.interior_surface = zeros(s_states,n_sur);%map surfaces to interior surface T state
ctf.surf_area = zeros(s_states,1);
ctf.subsurf_names = cell(s_states,1);
ctf.interior_absorb = zeros(n_sur,1);

n_ext = length(f_index({'Outdoors';'Ground'},surfaces.boundary));
ctf.exterior_surface = zeros(s_states,n_ext);%map surfaces to exterior surface T state
ctf.exterior_area = zeros(n_ext,1);
ctf.exterior_absorb = zeros(n_ext,1);
ctf.s_height = zeros(n_ext,1);
ctf.ground_cond = zeros(n_ext,1);
A_diag = zeros(s_states,3);% matrix for implicit calculation of heat conduction through surfaces (net flux from both sides = change in stored energy)
i = 0;
ext = 0;
for s = 1:1:n_sur
    j = f_index(surfaces.zone{s},z_names);%zone that it is touching
    ctf.map_sur2zone(j,s) = 1;%interior zone that this surface has convective heat transfer with
    C = surfaces.capacitance{s}*surfaces.area(s);
    ctf.interior_absorb(s) = surfaces.absorptance.thermal(s,1);
    i_nodes = length(surfaces.capacitance{s});
    K = surfaces.thermal_conductivity{s};
    if strcmp(surfaces.boundary{s},'Ground')
        K = K*6;
    end
    if first_side(s)
        ctf.interior_surface(i+1,s) = 1;
        for j = 1:i_nodes
            ctf.subsurf_names(i+j) = {strcat(surfaces.name{s},'_',num2str(j))};
        end
        A_diag(i+1:i+i_nodes,:) = implicit_matrices(K,surfaces.node_length{s});
        ctf.sur_state1(s) = i+1;
        ctf.capacitance(i+1:i+i_nodes) = C;
        ctf.surf_area(i+1:i+i_nodes) = surfaces.area(s);
        i = i+i_nodes;
    end
    switch surfaces.boundary{s}
        case {'Outdoors';'Ground'}
            ext = ext+1;
            if strcmp(surfaces.boundary{s},'Ground')
                ctf.ground_cond(ext,1) = K(i_nodes-1);
            end
            ctf.sur_state2(s) = i;
            ctf.exterior_surface(i,ext) = 1;
            ctf.exterior_area(ext) = surfaces.area(s);
            ctf.exterior_absorb(ext) = surfaces.absorptance.solar(s,2);
            ctf.s_height(ext) = mean(surfaces.vertices{s}(:,3));
        case 'Surface'
            b = f_index(surfaces.object{s},surfaces.name);
            if b>s%move second side up list to here so that A_diag matrice is collected together
                ctf.sur_state1(b) = i;
                ctf.subsurf_names(i) = surfaces.name(b);
                ctf.interior_surface(i,b) = 1;
            end
    end
end
A_diag(:,1) = [A_diag(2:end,1);0];%offset so sparse diag works correctly
A_diag(:,3) = [0;A_diag(1:end-1,3)];%offset so sparse diag works correctly
ctf.A = spdiags(A_diag,-1:1,s_states,s_states);
ctf.sur_state2 = nonzeros(ctf.sur_state2);
maxC = max(ctf.capacitance);
if min(ctf.capacitance)/max(ctf.capacitance)<1e-4
    ctf.capacitance(ctf.capacitance<5e-4*maxC) = 5e-4*maxC;%increase capacitance of air wall and the like to help with numerical solving
end
end%Ends function build_ctf

function A = implicit_matrices(K,dX)
i = length(K) + 1;
A = zeros(i,3);
for m = 1:i-1%conduction in both directions
    A(m,3) = -K(m)/dX(m);
    A(m+1,1) = -K(m)/dX(m);
end
A(:,2) = -sum(A,2);
end