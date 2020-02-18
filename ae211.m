%% Source & Vortex Panel Method for Batman Knife Shape 
% (https://www.throwninjastar.com/wp-content/uploads/2017/04/batman-throwing-star.jpg)
%Code by: Dr. Bilal A. Siddiqui(Asst. Prof., Mechanical Engineering, DHA Suffa University) & Deepak Prem
% MATLAB 2019a
%% Edited by: Apurva Nandan, Parag Sonkusare, Devendra Kharolia, Rahul Baudhh, Pranay Agrawal 
%% Group III

%% Batman Shape Definition
clc; clear all; close all;

Vinf=100;
n=500;%number of panels (n>16)
alpha=0;%angle of attack;
alpha = alpha*pi/180;
shape = 'batman';%'batman' or 'cylind' 

if(shape=='batman') 
    syms y(x)
    y(x) = piecewise((-17.54/0.9)<=x & x<=-8 , -0.5 *(.9*x + 17.54).^2, -8<x & x<=0, -0.6*(1.1*x + 7.5).^2 - 52.4438, 0<x & x<=8,-0.6*(-1.1*x + 7.5).^2 - 52.4438, (8<x & x<=17.54/0.9), -0.5*(-.9*x + 17.54).^2);

    dx=-2*17.54/(0.9*(n-14));
    X=17.54/0.9:dx:-17.54/0.9;
    Y=[];

    for itr=1:length(X)
        Y(itr)=y(X(itr));
    end

    X=[0,1,2,2.1,5,7,8,X,-8, -7, -5, -2.1, -2, -1, 0]';
    Y=[-9,-10,5,-68.5/3,-18,-11,-0.1,Y,-0.1,-11,-18,-68.5/3, 5,-10,-9]';
else
    R=10;
    dtheta=2*pi/n;
    theta=pi+pi/n:-dtheta:-pi-pi/n;%central angle
    X=R*cos(theta);
    Y=R*sin(theta);
end

%% Calculation of control points and other geometric parameters
for index=1:n
    %angle of flow with tangent of panel
    phi(index)=atan2((Y(index+1)-Y(index)),(X(index+1)-X(index)));
    %angle of flow with normal of panel
    beta(index)=phi(index)+pi/2;
    RHS(index)=sin(phi(index));
    midpoint_x(index)=(X(index+1)+X(index))/2;
    midpoint_y(index)=(Y(index+1)+Y(index))/2;
    %length of panel
    S(index)=sqrt((Y(index+1)-Y(index))^2+(X(index+1)-X(index))^2);
end

%% The Source Panel Method
for p=1:n
    xi=midpoint_x(p);
    yi=midpoint_y(p);
    for index=1:n
        if index == p
            continue
        end    
        Xj=X(index);
        Yj=Y(index);
        A=-(xi-Xj)*cos(phi(index))-(yi-Yj)*sin(phi(index));
        B=(xi-Xj)^2+(yi-Yj)^2;
        C=sin(phi(p)-phi(index));
        D=(yi-Yj)*cos(phi(p))-(xi-Xj)*sin(phi(p));
        E=sqrt(B-A^2);
        Sj=S(index);
        I(p,index)=C/2*log((Sj^2+2*A*Sj+B)/B)+(D-A*C)/E*(atan2((Sj+A),E)-atan2(A,E));
        J(p,index)=(D-A*C)/2/E*log((Sj^2+2*A*Sj+B)/B)-C*(atan2((Sj+A),E)-atan2(A,E));
    end
    F(p,1)=Vinf*cos(beta(p));
end
M=(I)/2/pi+eye(n)/2;

lambda=-inv(M)*F;

%Recoving velocity at the nodes
V=Vinf*sin(beta)+lambda'/2/pi*(J)';
Cp=1-(V/Vinf).^2;

fprintf('The sum of all sources by Source Panel Method is %f \n', (S*lambda));

%% The Vortex Panel Method
for i = 1:n
    for j = 1:n
        if (i == j)
            CN1(i,j) = -1;
            CN2(i,j) = 1 ;
            CT1(i,j) = 0.5*pi;
            CT2(i,j) = 0.5*pi;
        else
            A2 = - (midpoint_x(i) - X(j))*(cos(phi(j))) - (midpoint_y(i) - Y(j))*(sin(phi(j)));
            B2 = (midpoint_x(i) - X(j))^2 + (midpoint_y(i) - Y(j))^2;
            C2 = sin(phi(i) - phi(j));
            D2 = cos(phi(i) - phi(j));
            E2 = (midpoint_x(i) - X(j))*sin(phi(j)) - (midpoint_y(i) - Y(j))*cos(phi(j));
            F2 = log(1 + ((S(j))^2 + (2*A2*S(j))) / B2);
            G2 = atan2((E2*S(j)) , (B2 + A2*S(j)));
            P2 = ((midpoint_x(i) - X(j)) * sin(phi(i) - 2*phi(j))) + ((midpoint_y(i) - Y(j)) * cos(phi(i) - 2*phi(j)));
            Q2 = ((midpoint_x(i) - X(j)) * cos(phi(i) - 2*phi(j))) - ((midpoint_y(i) - Y(j)) * sin(phi(i) - 2*phi(j)));
            
            CN2(i,j) = D2 + ((0.5*Q2*F2)/S(j)) - ((A2*C2 + D2*E2)*(G2/S(j)));
            CN1(i,j) = 0.5*D2*F2 + C2*G2 - CN2(i,j);
            CT2(i,j) = C2 + ((0.5*P2*F2)/S(j)) + ((A2*D2 - C2*E2)*(G2/S(j)));
            CT1(i,j) = 0.5*C2*F2 - D2*G2 - CT2(i,j);
        end
    end
end
% Computation of Influence Coefficients
for i = 1:n
    AN(i,1) = CN1(i,1);
    AN(i,n+1) = CN2(1,n);
    AT(i,1) = CT1(i,1);
    AT(i,n+1) = CT2(i,n);
    for j = 2:n
        AN(i,j) = CN1(i,j) + CN2(i,j-1);
        AT(i,j) = CT1(i,j) + CT2(i,j-1);
    end
end
AN(n+1,1) = 1;
AN(n+1,n+1) = 1;
for j = 2:n
    AN(n+1,j) = 0;
end
RHS(n+1) = 0;

% Solve for Gamma and velocity/pressure
Gama = AN\(RHS');                % Solving for a syetem of linear equations
for i = 1:n
    V2(i) = cos(phi(i));
    for j = 1:n+1
        V2(i) = V2(i) + AT(i,j)*Gama(j);
        CP(i) = 1 - (V2(i))^2;
    end
end
% CP = [flip(CP(1:((n)/2))) CP((((n)/2)+1):end)];
circulation = sum(V2.*S);
fprintf('The sum of all sources by Vortex Panel Method is %f \n', sum(Gama(1:n)'.*(S)));
fprintf('The circulation associated with the shape by Vortex Panel Method is %e \n',circulation);

% Calculation of Lift Coefficient
fprintf('The lift coefficient of the shape evaluates to be %e by Vortex Panel Method\n',circulation/17.54*0.9);

% Plotting all the data
syms y1_graph(x);
syms y2_graph(x);
y1_graph(x) = piecewise((-17.54/0.9)<=x & x<-8 , 0.1/(-17.54/0.9+8)*(x+8)-0.1, -8<=x & x<-7, -11*x-88, -7<=x & x<-5,-3.5*(x+7)-11, -5<=x & x<-2.1, -5/3*(x+5)-18,-2.1<=x & x <-2,835/3*x-68.5/3+584.5 , -2<=x & x<-1,-15 *(x+1) - 10,-1<=x & x<0,-(-x - 1) - 10, (17.54/0.9)>=x & x>8 , 0.1/(17.54/0.9-8)*(x-8)-0.1, 8>=x & x>7, 11*x-88, 7>=x & x>5,-3.5*(-x+7)-11, 5>=x & x>2.1, -5/3*(-x+5)-18,2<x & x <=2.1,-835/3*x-68.5/3+584.5,2>x & x>1,-15 *(-x+1)-10,1>=x & x>=0,-(x - 1) - 10);
y2_graph(x) = piecewise((-17.54/0.9)<=x & x<=-8 , -0.5 *(.9*x + 17.54).^2, -8<x & x<0, -0.6*(1.1*x + 7.5).^2 - 52.4438, 0<=x & x<8,-0.6*(-1.1*x + 7.5).^2 - 52.4438, (8<=x & x<=17.54/0.9), -0.5*(-.9*x + 17.54).^2);
figure(1)
hold on
if(shape=='batman')
    h(1) = fplot(y2_graph(x), [-17.54/0.9,17.54/0.9], 'black','DisplayName','Exact Shape');
    h(2) = fplot(y1_graph(x), [-17.54/0.9,17.54/0.9], 'black','DisplayName','Exact Shape');
    h(3) = plot(X,Y,'r','DisplayName','Panel approximation');
    h(4) = plot(midpoint_x,midpoint_y,'g^','DisplayName','Control Points');
else
    ang=0:0.01:2*pi; 
    cir_x=R*cos(ang);
    cir_y=R*sin(ang);
    h(2) = plot(cir_y, cir_x, 'black','DisplayName','Exact Shape');
    h(3) = plot(X,Y,'r','DisplayName','Panel approximation');
    h(4) = plot(midpoint_x,midpoint_y,'g^','DisplayName','Control Points');
end
grid on;
grid minor;
legend(h(2:4));
figure(2)
plot(midpoint_x,Cp,'r');grid;
%legend('Cp (Source Panel Method)');
title('Cp Variation over the surface by Source Panel Method');
xlabel('X-coordinate');
ylabel('Coefficient of Pressure');
figure(3);
plot(midpoint_x,CP,'b');
% set(gca,'Ydir');
title('Cp Variation over the surface by Vortex Panel Method');
xlabel('X-coordinate');
ylabel('Coefficient of Pressure');
grid on;
grid minor;
hold off;
% close all;
