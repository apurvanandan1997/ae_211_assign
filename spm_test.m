%% Numerical source panel method:
clc;clear;close all;

%% Freestream flow:

q_inf = 1.0;
alpha = 0;
alpha = alpha*(pi/180);
u_inf = q_inf*cos(alpha);
v_inf = q_inf*sin(alpha);
%% Read Body coordinates:

filename = 'circular_cylinder.txt';
load(filename);

x_b = circular_cylinder(:,1);
y_b = circular_cylinder(:,2);
theta_b = circular_cylinder(:,3);
%theta_b(end) = '';
theta_b_deg = theta_b*(180/pi);
x_con = zeros(size(x_b)-1);
y_con = zeros(size(y_b)-1);


for i=1:1:length(x_b)-1
    
    x_con(i) = (x_b(i)+x_b(i+1))/2;
    y_con(i) = (y_b(i)+y_b(i+1))/2;
    theta_con(i) = (theta_b_deg(i)+theta_b_deg(i+1))/2; 
    theta_conrad(i) = (theta_b(i)+theta_b(i+1))/2; 
    
end

figure(1);hold on;
plot(x_b,y_b,'--bx')
plot(x_con,y_con,'ro')

%% RHS Vector:

radius = zeros(1,length(x_b)-1);
ang_cos = zeros(1,length(x_b)-1);
ang_sin = zeros(1,length(x_b)-1);
Normal_x = zeros(1,length(x_b)-1);
Normal_y = zeros(1,length(x_b)-1);
tangent_x = zeros(1,length(x_b)-1);
tangent_y = zeros(1,length(x_b)-1);
rhs_vec = zeros(1,length(x_b)-1);

for i=1:1:length(x_b)-1
    
    radius(i) = sqrt((x_b(i+1)-x_b(i))^2 + (y_b(i+1)-y_b(i))^2 );
    ang_cos(i) = (x_b(i+1)-x_b(i))/radius(i);
    ang_sin(i) = (y_b(i+1)-y_b(i))/radius(i);
    
    Normal_x(i) =  -ang_sin(i);
    Normal_y(i) = ang_cos(i);
    tangent_x(i) = ang_cos(i);
    tangent_y(i) = ang_sin(i);
    
    
    rhs_vec(i) = (Normal_x(i)*u_inf) + (Normal_y(i)*v_inf);
    
end

rhs_vec = -rhs_vec';

figure(2);hold on;
plot(x_b,y_b,'--b+');
quiver(x_con,y_con,Normal_x,Normal_y,0.1);
quiver(x_con,y_con,tangent_x,tangent_y,0.1)
%% LHS Matrix:

%% Preallocation phase:

radi = zeros(1,length(x_con));
cos_phii = zeros(size(radi));
sin_phii = zeros(size(radi));

radj = zeros(length(x_con),length(x_b)-1);
cos_phij = zeros(size(radj));
sin_phij = zeros(size(radj));
a = zeros(size(radj));
b = zeros(size(radj));
c = zeros(size(radj));
d = zeros(size(radj));
Sj = zeros(size(radj));
e = zeros(size(radj));
I1_ij = zeros(size(radj));
I2_ij = zeros(size(radj));
I3_ij = zeros(size(radj));
I_ij = zeros(size(radj));

cvind= zeros(size(radj));
dvind = zeros(size(radj));
I1ind_ij = zeros(size(radj));
I2ind_ij = zeros(size(radj));
I3ind_ij = zeros(size(radj));
Iind_ij = zeros(size(radj));


%% Compute phase:

for i  = 1:1:length(x_con)
    
    radi(i) = sqrt(((x_b(i+1)-x_b(i))^2)+((y_b(i+1)-y_b(i))^2));
    cos_phii(i) = (x_b(i+1)-x_b(i))/(radi(i));
    sin_phii(i) = (y_b(i+1)-y_b(i))/(radi(i));
    
    for j = 1:1:length(x_b)-1

        if(i ~= j)
        
      radj(i,j) = sqrt(((x_b(j+1)-x_b(j))^2)+((y_b(j+1)-y_b(j))^2));
      cos_phij(i,j) = (x_b(j+1)-x_b(j))/(radj(i,j));
      sin_phij(i,j) = (y_b(j+1)-y_b(j))/(radj(i,j));  
        
      a(i,j) = (-(x_con(i)-x_b(j))*cos_phij(i,j))-((y_con(i)-y_b(j))*sin_phij(i,j));
      b(i,j) = (x_con(i)-x_b(j))^2 + (y_con(i)-y_b(j))^2 ;
      c(i,j) = (sin_phii(i)*cos_phij(i,j))-(cos_phii(i)*sin_phij(i,j));
          
      d(i,j) = ((y_con(i)-y_b(j))*cos_phii(i))-((x_con(i)-x_b(j))*sin_phii(i)) ;
      Sj(i,j) = radj(i,j);
      e(i,j)  = ((x_con(i)-x_b(j))*(sin_phij(i,j)))-((y_con(i)-y_b(j))*(cos_phij(i,j)));
      
      I1_ij(i,j) = (c(i,j)/2)*log((Sj(i,j)*Sj(i,j)+(2*a(i,j)*Sj(i,j))+b(i,j))/(b(i,j)));
      I2_ij(i,j) = ((d(i,j)-(a(i,j)*c(i,j)))/(e(i,j)))*(atan((Sj(i,j)+a(i,j))/e(i,j)));
      I3_ij(i,j) = ((d(i,j)-(a(i,j)*c(i,j)))/(e(i,j)))*(atan(a(i,j)/e(i,j)));
      
      I_ij(i,j) = I1_ij(i,j)+I2_ij(i,j)-I3_ij(i,j);
      I_ij(i,j) = I_ij(i,j)*(1/(2*pi));
        
        else
            
            I_ij(i,j) = (1/2);
            
        end
      
        end
    
        
    
end

lambda = I_ij\rhs_vec;

figure(4);hold on;
plot(theta_con,lambda,':rx')

%% Induced velocity at panels:


for i = 1:1:length(x_con)
    
    for j = 1:1:length(x_con)
    
        if(i~=j)
        
           cvind(i,j) = -(sin_phii(i)*sin_phij(i,j)) - (cos_phii(i)*cos_phij(i,j)) ;
           dvind(i,j) = ((y_con(i)-y_b(j))*sin_phii(i))+((x_con(i)-x_b(j))*cos_phii(i)) ;
           
           I1ind_ij(i,j) = (cvind(i,j)/2)*log((Sj(i,j)*Sj(i,j)+(2*a(i,j)*Sj(i,j))+b(i,j))/(b(i,j)));
           I2ind_ij(i,j) = ((dvind(i,j)-(a(i,j)*cvind(i,j)))/(e(i,j)))*(atan((Sj(i,j)+a(i,j))/e(i,j)));
           I3ind_ij(i,j) = ((dvind(i,j)-(a(i,j)*cvind(i,j)))/(e(i,j)))*(atan(a(i,j)/e(i,j)));
      
           Iind_ij(i,j) = I1ind_ij(i,j)+I2ind_ij(i,j)-I3ind_ij(i,j);
           Iind_ij(i,j) = Iind_ij(i,j)*((lambda(j))/(2*pi));
            
        end
        
        if(i == j)
        
            
            Iind_ij(i,j) = (u_inf*tangent_x(i)) + (v_inf*tangent_y(i));
            
        end
    end
end

Vind = sum(Iind_ij,2);
Cp = 1-(Vind.^2)/(q_inf*q_inf);

Cpt = 1-4*(sin(theta_conrad).^2);
figure(5);hold on
plot(theta_con,Cp,'rx');
plot(theta_con,Cpt,'--b');

%% compute velocity field and streamlines:

%% grid creation:

delta_offset = 0.01;
Npt_rad = 70;
Npt_tan = 70;

radius = (x_b.^2)+(y_b.^2);
radius = max(max(radius));
rad_g = linspace(radius+delta_offset,3*radius,Npt_rad);
theta_g = linspace(0,2*pi,Npt_tan);
[rad_grid,theta_grid] = meshgrid(rad_g,theta_g);
x_grid = rad_grid.*cos(theta_grid);
y_grid = rad_grid.*sin(theta_grid);
figure(6);hold on;
plot(x_grid,y_grid,'rx');
hold on;
plot(x_b,y_b,'b+');

%% Compute velocity field:
Vx = zeros(size(x_grid));
Vy = zeros(size(x_grid));

for i = 1:1:size(x_grid,1)
    
    for j = 1:1:size(x_grid,2)
    
        Vx_t = 0;
        Vy_t = 0;
        
        for k=1:1:length(x_con)
        
            len_t = sqrt(((x_b(k+1)-x_b(k))^2)+((y_b(k+1)-y_b(k))^2));
            cphij = (x_b(k+1)-x_b(k))/(len_t);
            sphij = (y_b(k+1)-y_b(k))/(len_t);  
            at = -(x_grid(i,j)-x_b(k))*cphij -(y_grid(i,j)-y_b(k))*sphij;
            bt = (x_grid(i,j)-x_b(k))^2 + (y_grid(i,j)-y_b(k))^2;
            ct = -cphij;
            cty = -sphij;
            dt = (x_grid(i,j)-x_b(k));
            dty = (y_grid(i,j)-y_b(k));
            et = (x_grid(i,j)-x_b(k))*sphij -(y_grid(i,j)-y_b(k))*cphij;
            
            
           I1t = (ct/2)*log((len_t*len_t+(2*at*len_t)+bt)/(bt));
           I1ty = (cty/2)*log((len_t*len_t+(2*at*len_t)+bt)/(bt));
           
           I2t = ((dt-(at*ct))/(et))*(atan((len_t+at)/et));
           I2ty = ((dty-(at*cty))/(et))*(atan((len_t+at)/et));
           
           I3t = ((dt-(at*ct))/(et))*(atan(at/et));
           I3ty = ((dty-(at*cty))/(et))*(atan(at/et));
      
           It = I1t+I2t-I3t;
           It = It*((lambda(k))/(2*pi));
           
           Ity = I1ty+I2ty-I3ty;
           Ity = Ity*((lambda(k))/(2*pi));
           
           Vx_t = Vx_t + It;
           Vy_t = Vy_t + Ity;
            
        end
        
        Vx(i,j) = Vx_t;
        Vy(i,j) = Vy_t;
        
    end
end

%% Add free stream velocity components to solution:
Vx = Vx + u_inf;
Vy = Vy + v_inf;

Vtot = sqrt((Vx.^2)+(Vy.^2));
%%

figure(7);hold on;
plot(x_b,y_b,'b');
quiver(x_grid,y_grid,Vx,Vy);

figure(8);hold on;
surf(x_grid,y_grid,Vtot);

%% Plot tangential velocity at body surface: 

%% Put all negative tangential velocities to positive as this is due to
%% orientation/convention of surface tangent vectors on the lower side of thecylinder

Vindtemp = abs(Vind);

%%
V_surftan = Vtot(:,1);
figure(9);hold on;
plot(theta_g,V_surftan,'--rx');
plot(theta_conrad,Vindtemp,':b+')