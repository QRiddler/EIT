clear all;
clc;
x0 = [1,1];
f1 = @(x1,x2) 0.667*x1-1.33*x1*x2;
f2 = @(x1,x2) x2*x1-x2;
t0 = 0;
t_f = 10;
h =0.0001;
%% EULER METHOD
[y,t] = euler(f1,f2,t0,t_f,x0,h);
k = 0;
t = t0:h:t_f;
for i=1:length((t))
    x1(i)=y(k+1);
    x2(i)=y(k+2);
    k=k+2;
end
subplot(2,1,1)
plot(t,x1,'LineWidth',2)
legend('x1 by Euler method')
xlabel('t','FontSize',16,'FontName','Arial','FontWeight','bold')
ylabel('x1(t)','FontSize',16,'FontName','Arial','FontWeight','bold')
hold on
subplot(2,1,2)
plot(t,x2,'LineWidth',2)
legend('x2 by Euler method')
xlabel('t','FontSize',16,'FontName','Arial','FontWeight','bold')
ylabel('x2(t)','FontSize',16,'FontName','Arial','FontWeight','bold')
%% EULER METHOD FUNCTION
function [y,t] = euler(F1,F2,t0,t_f,y0,h);
y(1) = y0(1);
y(2) = y0(2);
t = t0:h:t_f;
k=0;
for i=1:length((t)-1)
    y(k+3)=y(k+1)+h*F1(y(k+1),y(k+2));
    y(k+4)=y(k+2)+h*F2(y(k+1),y(k+2));
    k=k+2;
end
end