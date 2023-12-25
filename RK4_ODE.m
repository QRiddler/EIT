clear all;
clc;
x0=[5,5,5]
f1=@(x1,x2,x3) 10*(x2-x1);
f2=@(x1,x2,x3) x1*(28-x3)-x2;
f3=@(x1,x2,x3) x1*x2-((8/3)*x3);
t0=0;
t_f=10;
h=0.01;
%% RK4
[y,t]=RK4ode(f1,f2,f3,t0,t_f,x0,h)
k=0;  
for i=1:length(t0:h:t_f)
    x1(i)=y(k+1);
    x2(i)=y(k+2);
    x3(i)=y(k+3);
    k = k+3;
end
subplot(3,1,1)
hold on
plot(t,x1,'linewidth',2)
legend('x1 by RK 4th order')
xlabel('t','Fontsize',16,'Fontname','Arial','fontweight','bold')
ylabel('x1(t)','Fontsize',16,'Fontname','Arial','fontweight','bold')
subplot(3,1,2)
plot(t,x2,'linewidth',2)
legend('x2 by RK 4th order')
xlabel('t','Fontsize',16,'Fontname','Arial','fontweight','bold')
ylabel('x2(t)','Fontsize',16,'Fontname','Arial','fontweight','bold')
subplot(3,1,3)
hold on
plot(t,x3,'linewidth',2)
legend('x3 by RK 4th order')
xlabel('t','Fontsize',16,'Fontname','Arial','fontweight','bold')
ylabel('x3(t)','Fontsize',16,'Fontname','Arial','fontweight','bold')
hold on
%% RK4 FUNCTION
function [y,t]=RK4ode(F1,F2,F3,t0,t_f,y0,h)
y(1)=y0(1);
y(2)=y0(2);
y(3)=y0(3);
t = t0:h:t_f;
j = 0;
for i=1:length((t)-1)
    k1(j+1)=F1(y(j+1),y(j+2),y(j+3));
    k1(j+2)=F1(y(j+1),y(j+2),y(j+3));
    k1(j+3)=F1(y(j+1),y(j+2),y(j+3));
    k2(j+1)=F1((y(j+1)+0.5*h*k1(j+1)),(y(j+2)+0.5*h*k1(j+2)),(y(j+3)+0.5*h*k1(j+3)));
    k2(j+2)=F2((y(j+1)+0.5*h*k1(j+1)),(y(j+2)+0.5*h*k1(j+2)),(y(j+3)+0.5*h*k1(j+3)));
    k2(j+3)=F3((y(j+1)+0.5*h*k1(j+1)),(y(j+2)+0.5*h*k1(j+2)),(y(j+3)+0.5*h*k1(j+3)));
    k3(j+1)=F1((y(j+1)+0.5*h*k2(j+1)),(y(j+2)+0.5*h*k2(j+2)),(y(j+3)+0.5*h*k2(j+3)));
    k3(j+2)=F2((y(j+1)+0.5*h*k2(j+1)),(y(j+2)+0.5*h*k2(j+2)),(y(j+3)+0.5*h*k2(j+3)));
    k3(j+3)=F3((y(j+1)+0.5*h*k2(j+1)),(y(j+2)+0.5*h*k2(j+2)),(y(j+3)+0.5*h*k2(j+3)));
    k4(j+1)=F1((y(j+1)+h*k3(j+1)),(y(j+2)+h*k3(j+2)),(y(j+3)+h*k2(j+3)));
    k4(j+2)=F2((y(j+1)+h*k3(j+1)),(y(j+2)+h*k3(j+2)),(y(j+3)+h*k2(j+3)));
    k4(j+3)=F3((y(j+1)+h*k3(j+1)),(y(j+2)+h*k3(j+2)),(y(j+3)+h*k2(j+3)));
    y(j+4)=y(j+1)+(h/6)*(k1(j+1)+2*k2(j+1)+2*k3(j+1)+k4(j+1));
    y(j+5)=y(j+2)+(h/6)*(k1(j+2)+2*k2(j+2)+2*k3(j+2)+k4(j+2));
    y(j+6)=y(j+3)+(h/6)*(k1(j+3)+2*k2(j+3)+2*k3(j+3)+k4(j+3));
    j=j+3;
end
end