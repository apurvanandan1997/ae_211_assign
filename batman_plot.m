syms y1(x);
syms y2(x);
y1(x) = piecewise((-17.54/0.9)<=x & x<-8 , 0.1/(-17.54/0.9+8)*(x+8)-0.1, -8<=x & x<-7, -11*x-88, -7<=x & x<-5,-3.5*(x+7)-11, -5<=x & x<-2.1, -5/3*(x+5)-18,-2.1<=x & x <-2,835/3*x-68.5/3+584.5 , -2<=x & x<-1,-15 *(x+1) - 10,-1<=x & x<0,-(-x - 1) - 10, (17.54/0.9)>=x & x>8 , 0.1/(17.54/0.9-8)*(x-8)-0.1, 8>=x & x>7, 11*x-88, 7>=x & x>5,-3.5*(-x+7)-11, 5>=x & x>2.1, -5/3*(-x+5)-18,2<x & x <=2.1,-835/3*x-68.5/3+584.5,2>x & x>1,-15 *(-x+1)-10,1>=x & x>=0,-(x - 1) - 10);
y2(x) = piecewise((-17.54/0.9)<=x & x<=-8 , -0.5 *(.9*x + 17.54).^2, -8<x & x<0, -0.6*(1.1*x + 7.5).^2 - 52.4438, 0<=x & x<8,-0.6*(-1.1*x + 7.5).^2 - 52.4438, (8<=x & x<=17.54/0.9), -0.5*(-.9*x + 17.54).^2);

hold on
fplot(y1(x), [-17.54/0.9,17.54/0.9], 'black')
fplot(y2(x), [-17.54/0.9,17.54/0.9], 'black')
hold off