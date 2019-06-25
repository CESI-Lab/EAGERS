function window_construct = window_glazing(window_construct,window_materials,mat,w)
switch window_construct.type{w}
    case 'SimpleGlazingSystem'
        window_construct = simple_glazing(window_construct,window_materials,mat,w);
    case 'Glazing'
        window_construct = glazing(window_construct,window_materials,mat,w);
end
end%Ends function glazing

function window = glazing(window,material,mat,w)
g = 0;
a = 0;
for i = 1:1:length(mat)
    j = mat(i);
    switch material.type{mat(i)}
        case 'Glazing'
            g = g+1;
            window.normal_transmittance(w,g) = material.solar_transmittance(j);
            window.emittance(w,2*g-1:2*g) = [material.emittance_front(j),material.emittance_back(j)];
            window.normal_reflectance(w,2*g-1:2*g) = [material.reflectance_front(j),material.reflectance_back(j)];
            window.thermal_conductivity(w,g) = material.thermal_conductivity(j);
            window.thickness(w,g) = material.thickness(j);
            if window.normal_transmittance(w,g)>0.645
                window.transmittance(w,1:5) = [0.0288, 1.460, -3.840, 3.355, -0.0015];
                window.reflectance(w,1:5) = [1.054, -2.532, 2.043, -0.563, 0.999];
            else
                window.transmittance(w,1:5) = [0.599, -0.05725, -2.341, 2.813, -0.002];
                window.reflectance(w,1:5) = [3.225, -7.862, 6.513, -1.868, 0.997];
            end
        case 'Gas'
            a = a+1;
            window.gap_thickness(w,a) = material.gap_thickness(j);
    end
end
end%Ends function glazing

function window = simple_glazing(window,window_materials,mat,w)
%conversion to single pane, see section 7.7, page 285 of reference manual
solar_heat_gain = window_materials.solar_heat_gain(mat(1));
u_factor = window_materials.u_factor(mat(1));
window.solar_heat_gain(w,1) = solar_heat_gain;
if isempty(window_materials.visible_transmittance(mat(1)))
    window.visible_transmittance(w,1) =nan;
else
    window.visible_transmittance(w,1) = window_materials.visible_transmittance(mat(1));
end
if u_factor<5.85
    window.interior_glazing_resistance_w(w,1) = 1/(0.359073*log(u_factor) + 6.949915);
else
    window.interior_glazing_resistance_w(w,1) = 1/(1.788041*u_factor - 2.886625);
end
window.exterior_glazing_resistance_w(w,1) = 1/(0.025342*u_factor + 29.163853);
window.bare_glass_resistance(w,1) = 1/u_factor - window.interior_glazing_resistance_w(w,1) - window.exterior_glazing_resistance_w(w,1);
if window.bare_glass_resistance(w,1)>7
    window.thickness(w,1) = 0.002;
else
    window.thickness(w,1) = 0.05914 - 0.00714/window.bare_glass_resistance(w,1);
end
window.thermal_conductivity(w,1) = window.thickness(w,1)/window.bare_glass_resistance(w,1);%%%

if solar_heat_gain<0.7206
    a =  0.939998*solar_heat_gain^2 + 0.20332*solar_heat_gain;
elseif solar_heat_gain>=0.7206
    a = 1.30415*solar_heat_gain - 0.30515;
end
if solar_heat_gain<=0.15
    b = 0.41040*solar_heat_gain;
elseif solar_heat_gain>0.15
    b = 0.085775*solar_heat_gain^2 + 0.963954*solar_heat_gain - 0.084958;
end
if u_factor>4.5
    window.normal_transmittance(w,1) = a;
elseif u_factor<3.4 
    window.normal_transmittance(w,1) = b;
else
    window.normal_transmittance(w,1) = (u_factor - 3.4)/1.1*(a-b)+b;
end

if isnan(window.visible_transmittance(w,1))
    window.visible_transmittance(w,1) = window.normal_transmittance(w,1);
end
SHGC_Tsol = solar_heat_gain - window.normal_transmittance(w,1);
a = 1/(29.436546*(SHGC_Tsol)^3 - 21.943415*(SHGC_Tsol)^2 + 9.945872*(SHGC_Tsol) + 7.426151);
b = 1/(199.8208128*(SHGC_Tsol)^3 - 90.639733*(SHGC_Tsol)^2 + 19.737055*(SHGC_Tsol) + 6.766575);
c = 1/(2.225824*(SHGC_Tsol) + 20.57708);
d = 1/(5.763355*(SHGC_Tsol) + 20.541528);
if u_factor>4.5
    window.interior_glazing_resistance_s(w,1) = a;
    window.exterior_glazing_resistance_s(w,1) = c;
elseif u_factor<3.4 
    window.interior_glazing_resistance_s(w,1) = b;
    window.exterior_glazing_resistance_s(w,1) = d;
else
    window.interior_glazing_resistance_s(w,1) = (u_factor - 3.4)/1.1*(a-b)+b;
    window.exterior_glazing_resistance_s(w,1) = (u_factor - 3.4)/1.1*(c-d)+d;
end
inward_fraction = (window.exterior_glazing_resistance_s(w,1) + 0.5*window.bare_glass_resistance(w,1))/(window.exterior_glazing_resistance_s(w,1) + window.bare_glass_resistance(w,1) + window.interior_glazing_resistance_s(w,1));
window.normal_reflectance(w,1:2) = 1 - window.normal_transmittance(w,1) - SHGC_Tsol/inward_fraction;
window.normal_reflectance(w,3:4) = 0;%2nd pane
window.emittance(w,1:2) = 0.84;%front side, back side
window.emittance(w,3:4) = 0;%2nd pane
window.long_wave_transmittance(w,1) = 0;
window.visible_reflectance(w,1) = -0.0622*window.visible_transmittance(w,1)^3 + 0.4277*window.visible_transmittance(w,1)^2 - 0.4169*window.visible_transmittance(w,1) + 0.2399;%front side
window.visible_reflectance(w,2) = -0.7409*window.visible_transmittance(w,1)^3 + 1.6531*window.visible_transmittance(w,1)^2 - 1.2299*window.visible_transmittance(w,1) + 0.4547;%baack side
window.visible_reflectance(w,3:4) = 0;%2nd pane
[window.transmittance(w,1:5),window.reflectance(w,1:5)] = transmit_reflect_coef(solar_heat_gain,u_factor);%%%%
window.gap_thickness(w,1) = 0;
window.normal_transmittance(w,2) = 0;%placeholder for 2nd pane
end%Ends function simple_glazing

function [transmittance,reflectance] = transmit_reflect_coef(SHGC,U)
%% Is not currently Used
%Details in engineering reference pg 287 and http://gaia.lbl.gov/btech/papers/2804.pdf
T = [1.470E-02 1.486E+00 -3.852E+00 3.355E+00 -1.474E-03;
     5.546E-01 3.563E-02 -2.416E+00 2.831E+00 -2.037E-03;
     7.709E-01 -6.383E-01 -1.576E+00 2.448E+00 -2.042E-03;
     3.462E-01 3.963E-01 -2.582E+00 2.845E+00 -2.804E-04;
     2.883E+00 -5.873E+00 2.489E+00 1.510E+00 -2.577E-03;
     3.025E+00 -6.366E+00 3.137E+00 1.213E+00 -1.367E-03;
     3.229E+00 -6.844E+00 3.535E+00 1.088E+00 -2.891E-03;
     3.334E+00 -7.131E+00 3.829E+00 9.766E-01 -2.952E-03;
     3.146E+00 -6.855E+00 3.931E+00 7.860E-01 -2.934E-03;
     3.744E+00 -8.836E+00 6.018E+00 8.407E-02 4.825E-04;];

R = [1.632E+01 -5.782E+01 7.924E+01 -5.008E+01 1.334E+01;
     4.048E+01 -1.193E+02 1.348E+02 -7.097E+01 1.611E+01;
     5.749E+01 -1.645E+02 1.780E+02 -8.875E+01 1.884E+01;
     5.714E+00 -1.667E+01 1.863E+01 -9.756E+00 3.074E+00;
     -5.488E-01 -6.498E+00 2.120E+01 -2.097E+01 7.814E+00;
     4.290E+00 -1.267E+01 1.466E+01 -8.153E+00 2.871E+00;
     2.174E+01 -6.444E+01 7.489E+01 -4.179E+01 1.062E+01;
     4.341E+00 -1.280E+01 1.478E+01 -8.203E+00 2.879E+00;
     4.136E+01 -1.178E+02 1.276E+02 -6.437E+01 1.426E+01;
     4.490E+00 -1.266E+01 1.397E+01 -7.501E+00 2.693E+00;]; 

T_fghi = .25*T(6,:) + .25*T(7,:) + .25*T(8,:) + .25*T(9,:); 
R_fghi = .25*R(6,:) + .25*R(7,:) + .25*R(8,:) + .25*R(9,:);
T_bdcd = .25*T(2,:) + .25*T(3,:) + .5*T(4,:);
R_bdcd = .25*R(2,:) + .25*R(3,:) + .5*R(4,:);
T_fh = 0.5*T(6,:) + 0.5*T(8,:);
R_fh = 0.5*R(6,:) + 0.5*R(8,:);
if SHGC>.65 
    if U>4.54
        transmittance = T(1,:);
        reflectance = R(1,:);
    elseif U>3.41
        a = (U-3.41)/1.13;
        transmittance = a*T(1,:) + (1-a)*T(5,:);    
        reflectance = a*R(1,:) + (1-a)*R(5,:);
    else%1,4,11
        transmittance = T(5,:);
        reflectance = R(5,:);       
    end
elseif SHGC>.6 
    b = (SHGC-.6)/.05;
    if U>4.54
        transmittance = b*T(1,:) + (1-b)*T_bdcd;  
        reflectance = b*R(1,:) + (1-b)*R_bdcd; 
    elseif U>3.41
        a = (U-3.41)/1.13;
        transmittance = a*(b*T(1,:) + (1-b)*T_bdcd) + (1-a)*T(5,:);    
        reflectance = a*(R(1,:) + (1-b)*R_bdcd) + (1-a)*R(5,:);
    else
        transmittance = T(5,:);
        reflectance = R(5,:);
    end
elseif SHGC>.55 
    if U>4.54
        transmittance = T_bdcd;  
        reflectance = R_bdcd; 
    elseif U>3.41
        a = (U-3.41)/1.13;
        transmittance = a*T_bdcd + (1-a)*T(5,:);    
        reflectance = a*R_bdcd + (1-a)*R(5,:);
    else
        transmittance = T(5,:);
        reflectance = R(5,:);
    end
elseif SHGC>.5 
    b = (SHGC-.5)/.05;
    if U>4.54
        transmittance = T_bdcd;  
        reflectance = R_bdcd; 
    elseif U>3.41
        a = (U-3.41)/1.13;
        transmittance = a*T_bdcd + (1-a)*(b*T(5,:) + (1-b)*T_fghi);    
        reflectance = a*R_bdcd + (1-a)*(b*R(5,:) + (1-b)*R_fghi);    
    elseif U>1.7
        transmittance = b*T(5,:) + (1-b)*T_fghi;
        reflectance = b*R(5,:) + (1-b)*R_fghi;
    elseif U>1.42
        a = (U-1.42)/0.28;
        transmittance = a*(b*T(5,:) + (1-b)*T_fghi) + (1-a)*T(5,:);    
        reflectance = a*(b*R(5,:) + (1-b)*R_fghi) + (1-a)*R(5,:);   
    else
        transmittance = T(5,:);
        reflectance = R(5,:);
    end
elseif SHGC>.45 
    if U>4.54
        transmittance = T_bdcd;  
        reflectance = R_bdcd; 
    elseif U>3.41
        a = (U-3.41)/1.13;
        transmittance = a*T_bdcd + (1-a)*T_fghi;    
        reflectance = a*R_bdcd + (1-a)*R_fghi;    
    elseif U>1.7
        transmittance = T_fghi;
        reflectance = R_fghi;
    elseif U>1.42
        a = (U-1.42)/0.28;
        transmittance = a*T_fghi + (1-a)*T(5,:);    
        reflectance = a*R_fghi + (1-a)*R(5,:);   
    else
        transmittance = T(5,:);
        reflectance = R(5,:);
    end
elseif SHGC>.35 
    if U>4.54
        b = (SHGC-.3)/.15;
        transmittance = b*T_bdcd + (1-b)*T(4,:);  
        reflectance = b*R_bdcd + (1-b)*R(4,:);  
    elseif U>3.41
        a = (U-3.41)/1.13;
        b = (SHGC-.3)/.15;
        transmittance = a*(b*T_bdcd + (1-b)*T(4,:)) + (1-a)*T_fghi;    
        reflectance = a*(b*R_bdcd + (1-b)*R(4,:)) + (1-a)*R_fghi;    
    elseif U>1.7
        transmittance = T_fghi;
        reflectance = R_fghi;
    elseif U>1.42
        b = (SHGC-.35)/.1;
        a = (U-1.42)/0.28;
        transmittance = a*T_fghi + (1-a)*(b*T(5,:) + (1-b)*T(10,:));    
        reflectance = a*R_fghi + (1-a)*(b*R(5,:) + (1-b)*R(10,:));   
    else
        b = (SHGC-.35)/.1;
        transmittance = b*T(5,:) + (1-b)*T(10,:);
        reflectance = b*R(5,:) + (1-b)*R(10,:);
    end
elseif SHGC>.3 
    if U>4.54
        b = (SHGC-.3)/.15;
        transmittance = b*T_bdcd + (1-b)*T(4,:);  
        reflectance = b*R_bdcd + (1-b)*R(4,:);  
    elseif U>3.41
        a = (U-3.41)/1.13;
        b = (SHGC-.3)/.15;
        transmittance = a*(b*T_bdcd + (1-b)*T(4,:)) + (1-a)*T_fghi;    
        reflectance = a*(b*R_bdcd + (1-b)*R(4,:)) + (1-a)*R_fghi;    
    elseif U>1.7
        transmittance = T_fghi;
        reflectance = R_fghi;
    elseif U>1.42
        a = (U-1.42)/0.28;
        transmittance = a*T_fghi + (1-a)*T(10,:);    
        reflectance = a*R_fghi + (1-a)*R(10,:);   
    else
        transmittance = T(10,:);
        reflectance = R(10,:);
    end    
elseif SHGC>.25 
    b = (SHGC-.25)/.05;
    if U>4.54
        transmittance = T(4,:);  
        reflectance = R(4,:);  
    elseif U>3.41
        a = (U-3.41)/1.13;
        transmittance = a*T(4,:) + (1-a)*(b*T_fghi + (1-b)*T_fh);    
        reflectance = a*R(4,:) + (1-a)*(b*R_fghi + (1-b)*R_fh);    
    elseif U>1.7
        transmittance = b*T_fghi + (1-b)*T_fh;
        reflectance = b*R_fghi + (1-b)*R_fh;
    elseif U>1.42
        a = (U-1.42)/0.28;
        transmittance = a*(T_fghi + (1-b)*T_fh) + (1-a)*T(10,:);    
        reflectance = a*(b*R_fghi + (1-b)*R_fh) + (1-a)*R(10,:);   
    else
        transmittance = T(10,:);
        reflectance = R(10,:);
    end  
else
    if U>4.54
        transmittance = T(4,:);  
        reflectance = R(4,:);  
    elseif U>3.41
        a = (U-3.41)/1.13;
        transmittance = a*T(4,:) + (1-a)*T_fh;    
        reflectance = a*R(4,:) + (1-a)*R_fh;    
    elseif U>1.7
        transmittance = T_fh;
        reflectance = R_fh;
    elseif U>1.42
        a = (U-1.42)/0.28;
        transmittance = a*(T_fh) + (1-a)*T(10,:);    
        reflectance = a*(R_fh) + (1-a)*R(10,:);   
    else
        transmittance = T(10,:);
        reflectance = R(10,:);
    end  
end
end%Ends funtion transmit_reflect_coef