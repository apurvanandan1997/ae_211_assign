function [ CL, Gamma, Cp, x, y] = Vortex_Panel( xb,yb,alpha )

%Define number of panels
M = length(xb)-1;

%Convert angle of attack from degrees to radians
alpha = alpha*pi/180;

%Pre-allocate arrays for panel location coordinates
x = zeros(1,M);
y = zeros(1,M);
s = zeros(1,M);
theta = zeros(1,M);
cosine = zeros(1,M);
sine = zeros(1,M);
RHS = zeros(1,M+1);
for i = 1:M
    ip1 = i+1;
    %Calculate the midpoint between the ith node and i+1th node
    x(i) = 0.5*(xb(i)+xb(ip1));
    y(i) = 0.5*(yb(i)+yb(ip1));
    s(i) = sqrt((xb(ip1)-xb(i))^2 +(yb(ip1)-yb(i))^2);
    theta(i) = atan2((yb(ip1)-yb(i)),(xb(ip1)-xb(i)));
    sine(i) = sin(theta(i));
    cosine(i) = cos(theta(i));
    RHS(i) = sin(theta(i)-alpha);
end

CN1 = zeros(M,M);
CN2 = zeros(M,M);
CT1 = zeros(M,M);
CT2 = zeros(M,M);
for i=1:M
   for j=1:M
      if i == j
         CN1(i,j) = -1.0;
         CN2(i,j) = 1.0;
         CT1(i,j) = 0.5*pi;
         CT2(i,j) = 0.5*pi;
      else
         A = -(x(i)-xb(j))*cosine(j) - (y(i)-yb(j))*sine(j);
         B = (x(i)-xb(j))^2 + (y(i)-yb(j))^2;
         C = sin(theta(i)-theta(j));
         D = cos(theta(i)-theta(j));
         E = (x(i)-xb(j))*sine(j) - (y(i)-yb(j))*cosine(j);
         F = log(1.0 + s(j)*(s(j)+2.*A)/B);
         G = atan2(E*s(j),B+A*s(j));
         P = (x(i)-xb(j))*sin(theta(i)-2*theta(j))+...
             (y(i)-yb(j))*cos(theta(i)-2*theta(j));
         Q = (x(i)-xb(j))*cos(theta(i)-2*theta(j))-...
             (y(i)-yb(j))*sin(theta(i)-2*theta(j));
         CN2(i,j) = D + 0.5*Q*F/s(j) - (A*C+D*E)*G/s(j);
         CN1(i,j) = 0.5*D*F + C*G - CN2(i,j);
         CT2(i,j) = C + 0.5*P*F/s(j) + (A*D-C*E)*G/s(j);
         CT1(i,j) = 0.5*C*F - D*G - CT2(i,j);
      end
   end
end

AN = zeros(M+1,M+1);
AT = zeros(M+1,M+1);
for i=1:M
    AN(i,1) = CN1(i,1);
    AN(i,M+1) = CN2(i,M);
    AT(i,1) = CT1(i,1);
    AT(i,M+1) = CT2(i,M);
    for j=2:M
       AN(i,j) = CN1(i,j) + CN2(i,j-1);
       AT(i,j) = CT1(i,j) + CT2(i,j-1);
    end
end
AN(M+1,1) = 1.0;
AN(M+1,M+1) = 1.0;

for j=2:M
   AN(M+1,j) = 0.0; 
end
RHS(M+1) = 0.0;

%Calculate vorticity
gamma = AN\RHS';
fprintf('The sum of all sources by Vortex Panel Method is %f \n', gamma(1:M)'.*(s));

%Calculate coefficient of pressure
V = zeros(1,M);
Cp = zeros(1,M);
for i=1:M
   V(i) = cos(theta(i) - alpha);
    for j=1:M+1
      V(i) = V(i) + AT(i,j)*gamma(j);
      Cp(i) = 1.0 - V(i)^2;
    end
end

%Calculate Coefficient of Lift
Gamma = 0;
for i=1:M
    Gamma = Gamma+V(i)*s(i);%/(max(xb)-min(xb));
end
CL = 2*Gamma/(max(xb)-min(xb));

%Plot Cp
plot(x/(max(x)-min(x)),Cp,'LineWidth',2);
set(gca ,'ydir','reverse');
xlabel('x/c');
ylabel('Coefficient of Pressure');
end
